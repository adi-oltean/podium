"""AST-based C99 emitter for static-subset kernels.

Supported subset (anything else raises EmitError — rejection IS the
spec): pure functions over float scalars and fixed-shape float arrays;
constant-index subscripts (1-D and 2-D); arithmetic with Python's
association order preserved (bit-exactness depends on it); `math.*`
calls from a whitelist (same libm as CPython on the host); `np.empty` /
`np.zeros` only as return-array allocation; if/else and conditional
expressions; calls to other emitted kernels (lowered through explicit
temporaries). No loops, no recursion, no heap, no exceptions.

Array-returning Python functions become void C functions with an `out`
parameter. Contracts from @contract render as ACSL `requires` clauses
and as the analyzer `[spec]` block documented in podium.verify.

Bit-exactness contract (tier-1 golden vectors): compile with
  gcc -std=c99 -O2 -ffp-contract=off
on x86-64 (SSE2 binary64). Expression trees are emitted with full
parenthesization in source order, so gcc evaluates exactly what CPython
evaluated; math.* maps to the same libm.
"""

from __future__ import annotations

import ast
import inspect
import textwrap
from dataclasses import dataclass, field
from typing import Callable

_MATH_FNS = {
    "sqrt", "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "fabs", "exp", "log", "hypot", "copysign", "floor", "fmod",
}


class EmitError(ValueError):
    """Source construct outside the static subset."""


@dataclass
class _FuncMeta:
    py_name: str
    c_name: str
    param_order: list[str] = field(default_factory=list)
    param_arrays: dict[str, int] = field(default_factory=dict)  # name->len
    out_shape: tuple[int, ...] | None = None  # None => returns scalar
    out_eye: bool = False  # identity-initialized return array (np.eye)
    contracts: dict = field(default_factory=dict)
    globals_: dict = field(default_factory=dict)


def _mangle(func: Callable) -> str:
    mod = func.__module__.rsplit(".", maxsplit=1)[-1]
    return f"podium_{mod}_{func.__name__}"


def _parse(func: Callable) -> ast.FunctionDef:
    src = textwrap.dedent(inspect.getsource(func))
    tree = ast.parse(src)
    node = tree.body[0]
    if not isinstance(node, ast.FunctionDef):
        raise EmitError(f"{func.__name__}: not a plain function")
    return node


def _const_index(node: ast.expr, fn: str) -> tuple[int, ...]:
    """Constant subscript index: int or tuple of ints."""
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return (node.value,)
    if isinstance(node, ast.Tuple):
        out = []
        for e in node.elts:
            if not (isinstance(e, ast.Constant) and isinstance(e.value, int)):
                raise EmitError(f"{fn}: non-constant subscript")
            out.append(e.value)
        return tuple(out)
    raise EmitError(f"{fn}: non-constant subscript")


def _analyze(func: Callable) -> tuple[_FuncMeta, ast.FunctionDef]:
    node = _parse(func)
    fn = func.__name__
    meta = _FuncMeta(py_name=fn, c_name=_mangle(func))
    meta.contracts = getattr(func, "__podium_contract__", {})
    meta.globals_ = getattr(func, "__wrapped__", func).__globals__
    meta.param_order = [a.arg for a in node.args.args]
    if node.args.vararg or node.args.kwarg or node.args.kwonlyargs \
            or node.args.defaults:
        raise EmitError(f"{fn}: only plain positional parameters supported")

    # array-ness + minimum lengths from constant subscripts (loop-var
    # indices contribute array-ness but no size bound)
    sizes: dict[str, int] = {}
    for sub in ast.walk(node):
        if isinstance(sub, ast.Subscript) and isinstance(sub.value, ast.Name):
            name = sub.value.id
            if name in meta.param_order:
                try:
                    idx = _const_index(sub.slice, fn)
                except EmitError:
                    sizes.setdefault(name, 0)
                    continue
                if len(idx) != 1:
                    raise EmitError(f"{fn}: 2-D parameter arrays unsupported")
                sizes[name] = max(sizes.get(name, 0), idx[0] + 1)
    meta.param_arrays = sizes

    # return shape from the allocation of the returned local
    ret = node.body[-1]
    if isinstance(ret, ast.Return) and isinstance(ret.value, ast.Name):
        ret_name = ret.value.id
        for st in node.body:
            if (isinstance(st, ast.Assign) and len(st.targets) == 1
                    and isinstance(st.targets[0], ast.Name)
                    and st.targets[0].id == ret_name
                    and isinstance(st.value, ast.Call)
                    and isinstance(st.value.func, ast.Attribute)
                    and st.value.func.attr in ("empty", "zeros", "eye")):
                arg = st.value.args[0]
                if st.value.func.attr == "eye":
                    if not isinstance(arg, ast.Constant):
                        raise EmitError(f"{fn}: non-constant eye size")
                    k = int(arg.value)
                    meta.out_shape = (k, k)
                    meta.out_eye = True
                elif isinstance(arg, ast.Constant):
                    meta.out_shape = (int(arg.value),)
                elif isinstance(arg, ast.Tuple):
                    meta.out_shape = tuple(
                        int(e.value) for e in arg.elts  # type: ignore[attr-defined]
                    )
                else:
                    raise EmitError(f"{fn}: non-constant array shape")
    return meta, node


class _Emitter(ast.NodeVisitor):
    def __init__(self, meta: _FuncMeta, registry: dict[str, _FuncMeta]):
        self.m = meta
        self.reg = registry
        self.lines: list[str] = []
        self.declared: set[str] = set()
        self.ret_name: str | None = None
        self.zero_init = False
        self.tmp_n = 0
        self.loop_vars: set[str] = set()

    # -- expressions ----------------------------------------------------
    def expr(self, e: ast.expr) -> str:
        if isinstance(e, ast.Constant):
            if isinstance(e.value, bool) or not isinstance(e.value, (int, float)):
                raise EmitError(f"{self.m.py_name}: unsupported constant {e.value!r}")
            return repr(float(e.value))
        if isinstance(e, ast.Name):
            if (e.id in self.declared or e.id in self.loop_vars
                    or e.id in self.m.param_order or e.id == self.ret_name):
                return e.id
            # module-level numeric constant (e.g. _TWO_PI): inline it
            g = self.m.globals_.get(e.id)
            if isinstance(g, (int, float)) and not isinstance(g, bool):
                return repr(float(g))
            raise EmitError(f"{self.m.py_name}: unknown name {e.id!r}")
        if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
            return f"(-{self.expr(e.operand)})"
        if isinstance(e, ast.BinOp):
            ops = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/"}
            for t, sym in ops.items():
                if isinstance(e.op, t):
                    return f"({self.expr(e.left)} {sym} {self.expr(e.right)})"
            raise EmitError(f"{self.m.py_name}: unsupported operator")
        if isinstance(e, ast.Subscript) and isinstance(e.value, ast.Name):
            name = ("out" if e.value.id == self.ret_name else e.value.id)
            return name + self._index_c(e.slice)
        if isinstance(e, ast.Compare):
            if len(e.ops) != 1:
                raise EmitError(f"{self.m.py_name}: chained comparison")
            sym = {ast.Gt: ">", ast.GtE: ">=", ast.Lt: "<", ast.LtE: "<=",
                   ast.Eq: "==", ast.NotEq: "!="}[type(e.ops[0])]
            return f"({self.expr(e.left)} {sym} {self.expr(e.comparators[0])})"
        if isinstance(e, ast.IfExp):
            return (f"({self.expr(e.test)} ? {self.expr(e.body)}"
                    f" : {self.expr(e.orelse)})")
        if isinstance(e, ast.Call):
            return self.call_expr(e)
        raise EmitError(f"{self.m.py_name}: unsupported expression "
                        f"{ast.dump(e)[:60]}")

    def _range_bound(self, st: ast.For) -> int:
        """Bounded `for _ in range(N)`: N a literal or a module-level
        integer constant — compile-time loop bounds, per the subset."""
        it = st.iter
        if not (isinstance(it, ast.Call) and isinstance(it.func, ast.Name)
                and it.func.id == "range" and len(it.args) == 1):
            raise EmitError(f"{self.m.py_name}: only range(N) loops")
        arg = it.args[0]
        if isinstance(arg, ast.Constant) and isinstance(arg.value, int):
            return arg.value
        if isinstance(arg, ast.Name):
            g = self.m.globals_.get(arg.id)
            if isinstance(g, int) and not isinstance(g, bool):
                return g
        raise EmitError(f"{self.m.py_name}: loop bound not compile-time")

    def _index_c(self, node: ast.expr) -> str:
        """Subscript index -> C brackets: constant ints, loop vars, or a
        tuple mixing them."""
        elts = node.elts if isinstance(node, ast.Tuple) else [node]
        parts = []
        for e in elts:
            if isinstance(e, ast.Constant) and isinstance(e.value, int):
                parts.append(f"[{e.value}]")
            elif isinstance(e, ast.Name) and e.id in self.loop_vars:
                parts.append(f"[{e.id}]")
            else:
                raise EmitError(f"{self.m.py_name}: unsupported subscript")
        return "".join(parts)

    def call_expr(self, e: ast.Call) -> str:
        # math.<fn>(...)
        if isinstance(e.func, ast.Attribute) and \
                isinstance(e.func.value, ast.Name) and e.func.value.id == "math":
            if e.func.attr not in _MATH_FNS:
                raise EmitError(f"{self.m.py_name}: math.{e.func.attr} not whitelisted")
            args = ", ".join(self.expr(a) for a in e.args)
            return f"{e.func.attr}({args})"
        # cross-kernel call returning an array -> lower to temp
        if isinstance(e.func, ast.Name) and e.func.id in self.reg:
            return self.lower_kernel_call(e)
        raise EmitError(f"{self.m.py_name}: unsupported call")

    def lower_kernel_call(self, e: ast.Call) -> str:
        callee = self.reg[e.func.id]  # type: ignore[union-attr]
        if callee.out_shape is None:
            args = ", ".join(self.expr(a) for a in e.args)
            return f"{callee.c_name}({args})"
        if len(callee.out_shape) != 1:
            raise EmitError(f"{self.m.py_name}: 2-D result in expression")
        tmp = f"_t{self.tmp_n}"
        self.tmp_n += 1
        self.lines.append(f"    double {tmp}[{callee.out_shape[0]}];")
        args = [self.expr(a) for a in e.args]
        self.lines.append(f"    {callee.c_name}({', '.join(args + [tmp])});")
        return tmp

    # -- statements -------------------------------------------------------
    def stmt(self, st: ast.stmt, indent: str = "    ") -> None:
        fn = self.m.py_name
        if isinstance(st, ast.Expr) and isinstance(st.value, ast.Constant):
            return  # docstring
        if isinstance(st, ast.Return):
            if isinstance(st.value, ast.Name) and st.value.id == self.ret_name:
                self.lines.append(f"{indent}return;")
                return
            if self.m.out_shape is None and st.value is not None:
                self.lines.append(f"{indent}return {self.expr(st.value)};")
                return
            raise EmitError(f"{fn}: unsupported return form")
        if isinstance(st, ast.Assign):
            if len(st.targets) != 1:
                raise EmitError(f"{fn}: multiple assignment")
            tgt = st.targets[0]
            # tuple unpack `a, b = e1, e2`: safe to sequence only when no
            # target name appears in any RHS element (checked)
            if isinstance(tgt, ast.Tuple):
                if not (isinstance(st.value, ast.Tuple)
                        and len(tgt.elts) == len(st.value.elts)
                        and all(isinstance(t, ast.Name) for t in tgt.elts)):
                    raise EmitError(f"{fn}: unsupported tuple assignment")
                names = {t.id for t in tgt.elts}  # type: ignore[union-attr]
                for rhs in st.value.elts:
                    for n_ in ast.walk(rhs):
                        if isinstance(n_, ast.Name) and n_.id in names:
                            raise EmitError(
                                f"{fn}: tuple assignment reads its targets")
                for t, rhs in zip(tgt.elts, st.value.elts):
                    self.stmt(ast.Assign(targets=[t], value=rhs), indent)
                return
            # the return-array allocation
            if isinstance(tgt, ast.Name) and isinstance(st.value, ast.Call) \
                    and isinstance(st.value.func, ast.Attribute) \
                    and st.value.func.attr in ("empty", "zeros", "eye") \
                    and isinstance(st.value.func.value, ast.Name) \
                    and st.value.func.value.id == "np":
                if self.ret_name is not None:
                    raise EmitError(f"{fn}: only one array allocation allowed")
                self.ret_name = tgt.id
                if st.value.func.attr in ("zeros", "eye"):
                    total = 1
                    for d in self.m.out_shape or ():
                        total *= d
                    self.lines.append(
                        f"{indent}for (int _i = 0; _i < {total}; _i++)"
                        f" ((double *)out)[_i] = 0.0;")
                if st.value.func.attr == "eye":
                    k = (self.m.out_shape or (0,))[0]
                    self.lines.append(
                        f"{indent}for (int _i = 0; _i < {k}; _i++)"
                        f" out[_i][_i] = 1.0;")
                return
            if isinstance(tgt, ast.Name):
                # array-valued kernel call assigned to a local: the local
                # IS the callee's out array
                if isinstance(st.value, ast.Call) \
                        and isinstance(st.value.func, ast.Name) \
                        and st.value.func.id in self.reg \
                        and self.reg[st.value.func.id].out_shape is not None:
                    callee = self.reg[st.value.func.id]
                    shape = callee.out_shape or ()
                    if len(shape) != 1:
                        raise EmitError(f"{fn}: 2-D call result unsupported")
                    if tgt.id in self.declared:
                        raise EmitError(f"{fn}: array local reassigned")
                    self.declared.add(tgt.id)
                    args = [self.expr(a) for a in st.value.args]
                    self.lines.append(
                        f"{indent}double {tgt.id}[{shape[0]}];")
                    self.lines.append(
                        f"{indent}{callee.c_name}"
                        f"({', '.join(args + [tgt.id])});")
                    return
                rhs = self.expr(st.value)
                if tgt.id in self.declared:
                    self.lines.append(f"{indent}{tgt.id} = {rhs};")
                else:
                    self.declared.add(tgt.id)
                    self.lines.append(f"{indent}double {tgt.id} = {rhs};")
                return
            if isinstance(tgt, ast.Subscript) and isinstance(tgt.value, ast.Name):
                name = "out" if tgt.value.id == self.ret_name else tgt.value.id
                lhs = name + self._index_c(tgt.slice)
                self.lines.append(f"{indent}{lhs} = {self.expr(st.value)};")
                return
            raise EmitError(f"{fn}: unsupported assignment target")
        if isinstance(st, ast.AugAssign):
            ops = {ast.Add: "+=", ast.Sub: "-=", ast.Mult: "*=",
                   ast.Div: "/="}
            sym = next((s for t, s in ops.items() if isinstance(st.op, t)),
                       None)
            if sym is None:
                raise EmitError(f"{fn}: unsupported augmented op")
            if isinstance(st.target, ast.Name):
                lhs = st.target.id
            elif isinstance(st.target, ast.Subscript) \
                    and isinstance(st.target.value, ast.Name):
                base = ("out" if st.target.value.id == self.ret_name
                        else st.target.value.id)
                lhs = base + self._index_c(st.target.slice)
            else:
                raise EmitError(f"{fn}: unsupported augmented target")
            self.lines.append(f"{indent}{lhs} {sym} {self.expr(st.value)};")
            return
        if isinstance(st, ast.For):
            bound = self._range_bound(st)
            var = st.target.id if isinstance(st.target, ast.Name) else None
            if var is None:
                raise EmitError(f"{fn}: loop target must be a name")
            self.loop_vars.add(var)
            self.lines.append(
                f"{indent}for (int {var} = 0; {var} < {bound}; {var}++) {{")
            for s in st.body:
                self.stmt(s, indent + "    ")
            self.lines.append(f"{indent}}}")
            self.loop_vars.discard(var)
            if st.orelse:
                raise EmitError(f"{fn}: for-else unsupported")
            return
        if isinstance(st, ast.If):
            self.lines.append(f"{indent}if {self.expr(st.test)} {{")
            for s in st.body:
                self.stmt(s, indent + "    ")
            if st.orelse:
                self.lines.append(f"{indent}}} else {{")
                for s in st.orelse:
                    self.stmt(s, indent + "    ")
            self.lines.append(f"{indent}}}")
            return
        raise EmitError(f"{fn}: unsupported statement {type(st).__name__}")


def _acsl(meta: _FuncMeta) -> str:
    """ACSL contract + analyzer [spec] block from @contract data."""
    if not meta.contracts:
        return ""
    req = []
    spec = []
    for name, iv in meta.contracts.items():
        if name in meta.param_arrays:
            n = meta.param_arrays[name]
            req.append(f"requires \\forall integer i; 0 <= i < {n} ==> "
                       f"{iv.lo!r} <= {name}[i] <= {iv.hi!r};")
        else:
            req.append(f"requires {iv.lo!r} <= {name} <= {iv.hi!r};")
        spec.append(f"[in, {iv.to_annotation()}] {name};")
    acsl = "/*@\n  " + "\n  ".join(req) + "\n*/\n"
    ann = "/* [spec] { " + " ".join(spec) + " } */\n"
    return acsl + ann


def _signature(meta: _FuncMeta) -> str:
    parts = []
    for p in meta.param_order:
        if p in meta.param_arrays:
            parts.append(f"const double {p}[{meta.param_arrays[p]}]")
        else:
            parts.append(f"double {p}")
    if meta.out_shape is not None:
        dims = "".join(f"[{d}]" for d in meta.out_shape)
        parts.append(f"double out{dims}")
        ret = "void"
    else:
        ret = "double"
    return f"{ret} {meta.c_name}({', '.join(parts) or 'void'})"


def emit_module(funcs: list[Callable], name: str = "podium_kernels") -> str:
    """Emit a self-contained C99 translation unit for the given kernels.

    Order matters only for readability; prototypes are emitted first so
    cross-kernel calls resolve regardless of definition order.
    """
    metas: dict[str, _FuncMeta] = {}
    nodes: dict[str, ast.FunctionDef] = {}
    for f in funcs:
        meta, node = _analyze(f)
        metas[meta.py_name] = meta
        nodes[meta.py_name] = node

    # propagate array-ness through cross-kernel calls: a parameter that
    # is never subscripted locally but is passed where a callee expects
    # an array IS an array of that size (two passes reach the fixpoint
    # for the call depths the subset allows)
    for _ in range(2):
        for py_name, meta in metas.items():
            for sub in ast.walk(nodes[py_name]):
                if not (isinstance(sub, ast.Call)
                        and isinstance(sub.func, ast.Name)
                        and sub.func.id in metas):
                    continue
                callee = metas[sub.func.id]
                for arg, pname in zip(sub.args, callee.param_order):
                    if isinstance(arg, ast.Name) \
                            and arg.id in meta.param_order \
                            and pname in callee.param_arrays:
                        size = callee.param_arrays[pname]
                        meta.param_arrays[arg.id] = max(
                            meta.param_arrays.get(arg.id, 0), size)

    out = [f"/* {name} — generated by podium.emit.cemit; DO NOT EDIT.",
           "   Compile: gcc -std=c99 -O2 -ffp-contract=off (SSE2/binary64);",
           "   bit-exactness vs CPython depends on those flags. */",
           "#include <math.h>", ""]
    for meta in metas.values():
        out.append(_signature(meta) + ";")
    out.append("")
    for py_name, meta in metas.items():
        em = _Emitter(meta, metas)
        for st in nodes[py_name].body:
            em.stmt(st)
        out.append(_acsl(meta) + _signature(meta) + " {")
        out.extend(em.lines)
        out.append("}")
        out.append("")
    return "\n".join(out)
