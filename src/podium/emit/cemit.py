"""AST-based C99 emitter for static-subset kernels.

Supported subset (anything else raises EmitError — rejection IS the
spec): pure functions over float scalars and fixed-shape float arrays;
constant-index subscripts (1-D and 2-D); arithmetic with Python's
association order preserved (bit-exactness depends on it); `math.*`
calls from a whitelist (same libm as CPython on the host); `np.empty` /
`np.zeros` only as return-array allocation; if/else and conditional
expressions; calls to other emitted kernels (lowered through explicit
temporaries); bounded compile-time `range(N)` loops. No recursion, no
heap, no exceptions.

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
    param_shapes: dict[str, tuple[int, ...]] = field(default_factory=dict)
    out_shape: tuple[int, ...] | None = None  # None => returns scalar
    out_eye: bool = False  # identity-initialized return array (np.eye)
    # tuple returns: ordered (local name, shape) -> out0, out1, ...
    outs: list[tuple[str, tuple[int, ...]]] = field(default_factory=list)
    allocs: dict[str, tuple[tuple[int, ...], str]] = field(
        default_factory=dict)  # local array name -> (shape, kind)
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

    # declared shapes (@shapes decorator) take precedence for params
    declared = getattr(func, "__podium_shapes__", {})
    for p, shp in declared.items():
        if p not in meta.param_order:
            raise EmitError(f"{fn}: @shapes names unknown parameter {p!r}")
        if len(shp) == 1:
            meta.param_arrays[p] = shp[0]
        meta.param_shapes[p] = tuple(shp)

    # array-ness + minimum lengths from constant subscripts (loop-var
    # indices contribute array-ness but no size bound)
    sizes: dict[str, int] = dict(meta.param_arrays)
    for sub in ast.walk(node):
        if isinstance(sub, ast.Subscript) and isinstance(sub.value, ast.Name):
            name = sub.value.id
            if name in meta.param_order and name not in meta.param_shapes:
                try:
                    idx = _const_index(sub.slice, fn)
                except EmitError:
                    sizes.setdefault(name, 0)
                    continue
                if len(idx) != 1:
                    raise EmitError(f"{fn}: 2-D parameter arrays need @shapes")
                sizes[name] = max(sizes.get(name, 0), idx[0] + 1)
    meta.param_arrays = sizes
    for p, ln in sizes.items():
        meta.param_shapes.setdefault(p, (ln,))

    # local array allocations (np.empty/zeros/eye)
    def _alloc_shape(call: ast.Call) -> tuple[tuple[int, ...], str]:
        kind = call.func.attr  # type: ignore[union-attr]
        arg = call.args[0]
        if kind == "eye":
            if not isinstance(arg, ast.Constant):
                raise EmitError(f"{fn}: non-constant eye size")
            return (int(arg.value), int(arg.value)), kind
        if isinstance(arg, ast.Constant):
            return (int(arg.value),), kind
        if isinstance(arg, ast.Tuple):
            return tuple(int(e.value) for e in arg.elts), kind  # type: ignore[attr-defined]
        raise EmitError(f"{fn}: non-constant array shape")

    for st in ast.walk(node):
        if (isinstance(st, ast.Assign) and len(st.targets) == 1
                and isinstance(st.targets[0], ast.Name)
                and isinstance(st.value, ast.Call)
                and isinstance(st.value.func, ast.Attribute)
                and isinstance(st.value.func.value, ast.Name)
                and st.value.func.value.id == "np"
                and st.value.func.attr in ("empty", "zeros", "eye")):
            meta.allocs[st.targets[0].id] = _alloc_shape(st.value)

    # shape-inference pass: arrays created by pure matrix EXPRESSIONS
    # (x_out = phi @ x) never appear as np allocations, but they can be
    # returned — infer their shapes statically so signatures resolve
    env = dict(meta.param_shapes)

    def _scan(body: list) -> None:
        for raw in body:
            st = raw
            if isinstance(raw, ast.AnnAssign) and raw.value is not None \
                    and isinstance(raw.target, ast.Name):
                st = ast.Assign(targets=[raw.target], value=raw.value)
            if isinstance(st, ast.Assign) and len(st.targets) == 1 \
                    and isinstance(st.targets[0], ast.Name):
                nm = st.targets[0].id
                if nm in meta.allocs:
                    env[nm] = meta.allocs[nm][0]
                    continue
                try:
                    shp = _ashape(st.value, env, meta.globals_)
                except EmitError:
                    shp = None
                if shp is not None:
                    env[nm] = shp
                    meta.allocs.setdefault(nm, (shp, "expr"))
            elif isinstance(st, (ast.For, ast.If)):
                _scan(st.body)
                _scan(list(getattr(st, "orelse", [])))

    _scan(node.body)

    # returns: single local array, tuple of locals, or scalar expression
    ret = node.body[-1]
    if isinstance(ret, ast.Return) and isinstance(ret.value, ast.Name) \
            and ret.value.id in meta.allocs:
        shape, kind = meta.allocs[ret.value.id]
        meta.out_shape = shape
        meta.out_eye = kind == "eye"
        meta.outs = [(ret.value.id, shape)]
    elif isinstance(ret, ast.Return) and isinstance(ret.value, ast.Tuple):
        for el in ret.value.elts:
            if not (isinstance(el, ast.Name) and el.id in meta.allocs):
                raise EmitError(f"{fn}: tuple return must name local arrays")
            meta.outs.append((el.id, meta.allocs[el.id][0]))
    return meta, node


def _ashape(e: ast.expr, env: dict, globals_: dict) -> tuple[int, ...] | None:
    """Static shape of an expression (None => scalar); raises on
    dimension mismatches so translation bugs die at emit time."""
    if isinstance(e, ast.Name):
        return env.get(e.id)
    if isinstance(e, ast.Attribute) and e.attr == "T":
        s = _ashape(e.value, env, globals_)
        if s is None:
            raise EmitError("transpose of a scalar")
        return tuple(reversed(s)) if len(s) == 2 else s
    if isinstance(e, ast.UnaryOp):
        return _ashape(e.operand, env, globals_)
    if isinstance(e, ast.BinOp):
        ls = _ashape(e.left, env, globals_)
        rs = _ashape(e.right, env, globals_)
        if isinstance(e.op, ast.MatMult):
            if ls is None or rs is None:
                raise EmitError("@ needs two arrays")
            if len(ls) == 2 and len(rs) == 1 and ls[1] == rs[0]:
                return (ls[0],)
            if len(ls) == 2 and len(rs) == 2 and ls[1] == rs[0]:
                return (ls[0], rs[1])
            raise EmitError(f"@ shape mismatch {ls} x {rs}")
        if ls is None and rs is None:
            return None
        if isinstance(e.op, (ast.Add, ast.Sub)):
            if ls is not None and rs is not None and ls == rs:
                return ls
            raise EmitError(f"elementwise +/- shape mismatch {ls} vs {rs}")
        if isinstance(e.op, (ast.Mult, ast.Div)):
            if ls is not None and rs is not None:
                if ls == rs and isinstance(e.op, ast.Mult):
                    return ls
                raise EmitError("array*array supported only same-shape")
            return ls if ls is not None else rs
        raise EmitError("unsupported array operator")
    return None


# libm name -> CORE-MATH correctly-rounded name, in cr mode
_CR_MAP = {"sin": "cr_sin", "cos": "cr_cos"}


class _Emitter(ast.NodeVisitor):
    def __init__(self, meta: _FuncMeta, registry: dict[str, _FuncMeta],
                 correctly_rounded: bool = False):
        self.m = meta
        self.reg = registry
        self.cr = correctly_rounded
        self.lines: list[str] = []
        self.declared: set[str] = set()
        self.env: dict[str, tuple[int, ...]] = dict(meta.param_shapes)
        # return arrays -> C out-parameter names
        if len(meta.outs) == 1:
            self.outs_map = {meta.outs[0][0]: "out"}
        else:
            self.outs_map = {n: f"out{i}" for i, (n, _s) in
                             enumerate(meta.outs)}
        self.tmp_n = 0
        self.idx_n = 0
        self.loop_vars: set[str] = set()
        self.loop_bounds: dict[str, int] = {}  # loop var -> exclusive bound

    def cn(self, name: str) -> str:
        return self.outs_map.get(name, name)

    def _dim_size(self, name: str, dim: int) -> int | None:
        """Declared size of dimension ``dim`` of array ``name``, or None if
        not known at emit time (parameters, locals, temporaries)."""
        shape: tuple[int, ...] | None = None
        if name in self.m.param_shapes:
            shape = self.m.param_shapes[name]
        elif name in self.m.allocs:
            shape = self.m.allocs[name][0]
        elif name in self.env:
            shape = self.env[name]
        if shape is None or dim >= len(shape):
            return None
        # a param indexed only by loop variables has no size inferred from
        # constant subscripts (recorded as 0); treat that as "unknown", not
        # as a real zero-length array.
        return shape[dim] or None

    def _fresh_idx(self) -> str:
        self.idx_n += 1
        return f"_i{self.idx_n}"

    def _int_typed(self, e: ast.expr) -> bool:
        """True if ``e`` emits as an int-typed C expression. Loop counters
        (declared ``int``) are the subset's only integers; ``/`` between two
        such operands is C integer division, which truncates and traps on a
        zero divisor rather than matching Python's float division."""
        if isinstance(e, ast.Name):
            return e.id in self.loop_vars
        if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
            return self._int_typed(e.operand)
        if isinstance(e, ast.BinOp) and isinstance(
                e.op, (ast.Add, ast.Sub, ast.Mult)):
            return self._int_typed(e.left) and self._int_typed(e.right)
        return False

    # -- expressions ----------------------------------------------------
    def expr(self, e: ast.expr) -> str:
        if isinstance(e, ast.Constant):
            if isinstance(e.value, bool) or not isinstance(e.value, (int, float)):
                raise EmitError(f"{self.m.py_name}: unsupported constant {e.value!r}")
            return repr(float(e.value))
        if isinstance(e, ast.Name):
            if (e.id in self.declared or e.id in self.loop_vars
                    or e.id in self.m.param_order or e.id in self.env
                    or e.id in self.outs_map):
                return self.cn(e.id)
            # module-level numeric constant (e.g. _TWO_PI): inline it
            g = self.m.globals_.get(e.id)
            if isinstance(g, (int, float)) and not isinstance(g, bool):
                return repr(float(g))
            raise EmitError(f"{self.m.py_name}: unknown name {e.id!r}")
        if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
            return f"(-{self.expr(e.operand)})"
        if isinstance(e, ast.BinOp):
            if isinstance(e.op, ast.Div) and self._int_typed(e.left) \
                    and self._int_typed(e.right):
                raise EmitError(
                    f"{self.m.py_name}: integer division of loop variables "
                    f"(C truncates toward zero and traps on a zero divisor; "
                    f"Python does true float division)")
            ops = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/"}
            for t, sym in ops.items():
                if isinstance(e.op, t):
                    return f"({self.expr(e.left)} {sym} {self.expr(e.right)})"
            raise EmitError(f"{self.m.py_name}: unsupported operator")
        if isinstance(e, ast.Subscript) and isinstance(e.value, ast.Name):
            return self.cn(e.value.id) + self._index_c(e.slice, e.value.id)
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

    def _index_c(self, node: ast.expr, arr: str | None = None) -> str:
        """Subscript index -> C brackets: constant ints, loop vars,
        loop-var + constant arithmetic, or tuples mixing them. Every index
        is checked to stay within the array's declared bounds where those
        are known, so the emitted C cannot read out of bounds where Python
        would raise IndexError."""
        elts = node.elts if isinstance(node, ast.Tuple) else [node]
        parts = []
        for dim, e in enumerate(elts):
            size = self._dim_size(arr, dim) if arr is not None else None
            if isinstance(e, ast.Constant) and isinstance(e.value, int):
                if size is not None and e.value >= size:
                    raise EmitError(
                        f"{self.m.py_name}: subscript [{e.value}] out of "
                        f"bounds for size {size}")
                parts.append(f"[{e.value}]")
            elif isinstance(e, ast.Name) and e.id in self.loop_vars:
                bound = self.loop_bounds.get(e.id)
                if size is not None and bound is not None and bound > size:
                    raise EmitError(
                        f"{self.m.py_name}: loop {e.id} in range({bound}) "
                        f"indexes past size {size}")
                parts.append(f"[{e.id}]")
            elif (isinstance(e, ast.BinOp)
                  and isinstance(e.op, (ast.Add, ast.Sub))
                  and isinstance(e.left, ast.Name)
                  and e.left.id in self.loop_vars
                  and isinstance(e.right, ast.Constant)
                  and isinstance(e.right.value, int)):
                # a loop variable ranges over [0, bound); a negative net
                # offset (e.g. i - 1) can make the index negative, which
                # Python wraps to the end but C reads out of bounds.
                off = (e.right.value if isinstance(e.op, ast.Add)
                       else -e.right.value)
                if off < 0:
                    raise EmitError(
                        f"{self.m.py_name}: subscript {e.left.id} - "
                        f"{e.right.value} may be negative (Python wraps to "
                        f"the end, C reads out of bounds)")
                # a positive offset can over-read; verify the maximum index
                # (bound-1)+off stays within the declared size. Relative
                # indexing requires both to be known so it can be proved safe.
                bound = self.loop_bounds.get(e.left.id)
                if bound is None or size is None:
                    raise EmitError(
                        f"{self.m.py_name}: offset subscript {e.left.id} + "
                        f"{off} needs a known loop bound and array size to be "
                        f"proved in bounds")
                if (bound - 1) + off >= size:
                    raise EmitError(
                        f"{self.m.py_name}: subscript {e.left.id} + {off} "
                        f"reaches index {(bound - 1) + off} >= size {size} "
                        f"(out of bounds in C)")
                parts.append(f"[{e.left.id}]" if off == 0
                             else f"[{e.left.id} + {off}]")
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
            fn = e.func.attr
            if self.cr and fn in _CR_MAP:
                fn = _CR_MAP[fn]
            return f"{fn}({args})"
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
        self.env[tmp] = callee.out_shape
        return tmp

    # -- array-expression lowering ---------------------------------------
    def _decl_array(self, shape: tuple[int, ...], indent: str) -> str:
        tmp = f"_a{self.tmp_n}"
        self.tmp_n += 1
        dims = "".join(f"[{d}]" for d in shape)
        self.lines.append(f"{indent}double {tmp}{dims};")
        self.env[tmp] = shape
        return tmp

    def _resolve(self, e: ast.expr, indent: str) -> tuple[str, tuple[int, ...], bool]:
        """Resolve an array-expression operand to (c_name, shape,
        transposed) — materializing sub-expressions as temporaries."""
        if isinstance(e, ast.Name) and e.id in self.env:
            return self.cn(e.id), self.env[e.id], False
        if isinstance(e, ast.Attribute) and e.attr == "T":
            base, shape, trans = self._resolve(e.value, indent)
            if len(shape) != 2:
                return base, shape, trans  # 1-D transpose is identity
            return base, shape, not trans
        # anything else: lower into a temp
        name = self.lower_array(e, None, indent)
        return name, self.env[name], False

    def _at(self, name: str, shape: tuple[int, ...], trans: bool,
            i: str, j: str = "") -> str:
        if len(shape) == 1:
            return f"{name}[{i}]"
        return f"{name}[{j}][{i}]" if trans else f"{name}[{i}][{j}]"

    def lower_array(self, e: ast.expr, dest: str | None, indent: str) -> str:
        """Emit code computing array expression `e` into `dest` (or a
        fresh temp); returns the destination C name. Loops are explicit
        with fixed bounds; summation is naive row-major order (NumPy's
        BLAS order differs, so matmul kernels are a relative-tolerance
        golden-vector class, documented in the tests)."""
        fn = self.m.py_name
        shape = _ashape(e, self.env, self.m.globals_)
        assert shape is not None     # _ashape is non-None for every array-shaped caller
        if isinstance(e, ast.Name):
            src = self.cn(e.id)
            if dest is None or dest == src:
                return src
            d = dest
            self._copy(src, d, shape, indent)
            return d
        if isinstance(e, ast.Attribute) and e.attr == "T" and len(shape) == 1:
            return self.lower_array(e.value, dest, indent)
        d = dest if dest is not None else self._decl_array(shape, indent)
        if isinstance(e, ast.BinOp) and isinstance(e.op, ast.MatMult):
            a, ash, at_ = self._resolve(e.left, indent)
            b, bsh, bt = self._resolve(e.right, indent)
            # effective shapes after transposition
            ea = tuple(reversed(ash)) if at_ else ash
            eb = tuple(reversed(bsh)) if bt else bsh
            i, j, k = self._fresh_idx(), self._fresh_idx(), self._fresh_idx()
            if len(eb) == 1:  # (m,k) @ (k,)
                self.lines.append(
                    f"{indent}for (int {i} = 0; {i} < {ea[0]}; {i}++) {{")
                self.lines.append(f"{indent}    double _acc = 0.0;")
                self.lines.append(
                    f"{indent}    for (int {k} = 0; {k} < {ea[1]}; {k}++)"
                    f" _acc += {self._at(a, ash, at_, i, k)} * {b}[{k}];")
                self.lines.append(f"{indent}    {d}[{i}] = _acc;")
                self.lines.append(f"{indent}}}")
            else:  # (m,k) @ (k,n)
                self.lines.append(
                    f"{indent}for (int {i} = 0; {i} < {ea[0]}; {i}++)"
                    f" for (int {j} = 0; {j} < {eb[1]}; {j}++) {{")
                self.lines.append(f"{indent}    double _acc = 0.0;")
                self.lines.append(
                    f"{indent}    for (int {k} = 0; {k} < {ea[1]}; {k}++)"
                    f" _acc += {self._at(a, ash, at_, i, k)}"
                    f" * {self._at(b, bsh, bt, k, j)};")
                self.lines.append(f"{indent}    {d}[{i}][{j}] = _acc;")
                self.lines.append(f"{indent}}}")
            return d
        if isinstance(e, ast.BinOp) and isinstance(e.op, (ast.Add, ast.Sub,
                                                          ast.Mult, ast.Div)):
            sym = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*",
                   ast.Div: "/"}[type(e.op)]
            ls = _ashape(e.left, self.env, self.m.globals_)
            rs = _ashape(e.right, self.env, self.m.globals_)
            i, j = self._fresh_idx(), self._fresh_idx()

            def elem(side: ast.expr,
                     sshape: tuple[int, ...] | None,
                     name_trans: tuple[str, tuple[int, ...], bool] | None,
                     ) -> str:
                if sshape is None or name_trans is None:
                    return self.expr(side)
                nm, sh, tr = name_trans
                return self._at(nm, sh, tr, i, j if len(shape) == 2 else "")

            lnt = self._resolve(e.left, indent) if ls is not None else None
            rnt = self._resolve(e.right, indent) if rs is not None else None
            if len(shape) == 1:
                self.lines.append(
                    f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)")
                self.lines.append(
                    f"{indent}    {d}[{i}] = {elem(e.left, ls, lnt)} {sym} "
                    f"{elem(e.right, rs, rnt)};")
            else:
                self.lines.append(
                    f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)"
                    f" for (int {j} = 0; {j} < {shape[1]}; {j}++)")
                self.lines.append(
                    f"{indent}    {d}[{i}][{j}] = {elem(e.left, ls, lnt)}"
                    f" {sym} {elem(e.right, rs, rnt)};")
            return d
        if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
            src, sh, tr = self._resolve(e.operand, indent)
            i, j = self._fresh_idx(), self._fresh_idx()
            if len(shape) == 1:
                self.lines.append(
                    f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)"
                    f" {d}[{i}] = -{self._at(src, sh, tr, i)};")
            else:
                self.lines.append(
                    f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)"
                    f" for (int {j} = 0; {j} < {shape[1]}; {j}++)"
                    f" {d}[{i}][{j}] = -{self._at(src, sh, tr, i, j)};")
            return d
        if isinstance(e, ast.Attribute) and e.attr == "T":
            src, sh, tr = self._resolve(e, indent)
            i, j = self._fresh_idx(), self._fresh_idx()
            self.lines.append(
                f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)"
                f" for (int {j} = 0; {j} < {shape[1]}; {j}++)"
                f" {d}[{i}][{j}] = {self._at(src, sh, tr, i, j)};")
            return d
        raise EmitError(f"{fn}: unsupported array expression "
                        f"{ast.dump(e)[:60]}")

    def _copy(self, src: str, dst: str, shape: tuple[int, ...],
              indent: str) -> None:
        i, j = self._fresh_idx(), self._fresh_idx()
        if len(shape) == 1:
            self.lines.append(
                f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)"
                f" {dst}[{i}] = {src}[{i}];")
        else:
            self.lines.append(
                f"{indent}for (int {i} = 0; {i} < {shape[0]}; {i}++)"
                f" for (int {j} = 0; {j} < {shape[1]}; {j}++)"
                f" {dst}[{i}][{j}] = {src}[{i}][{j}];")

    # -- statements -------------------------------------------------------
    def stmt(self, st: ast.stmt, indent: str = "    ") -> None:
        fn = self.m.py_name
        if isinstance(st, ast.Expr) and isinstance(st.value, ast.Constant):
            return  # docstring
        if isinstance(st, ast.AnnAssign):
            if st.value is None:
                raise EmitError(f"{fn}: bare annotation unsupported")
            st = ast.Assign(targets=[st.target], value=st.value)
        if isinstance(st, ast.Return):
            if isinstance(st.value, ast.Name) and st.value.id in self.outs_map:
                self.lines.append(f"{indent}return;")
                return
            if isinstance(st.value, ast.Tuple) and all(
                    isinstance(el, ast.Name) and el.id in self.outs_map
                    for el in st.value.elts):
                self.lines.append(f"{indent}return;")
                return
            if not self.m.outs and st.value is not None:
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
            # explicit np allocation of a local/return array
            if isinstance(tgt, ast.Name) and isinstance(st.value, ast.Call) \
                    and isinstance(st.value.func, ast.Attribute) \
                    and st.value.func.attr in ("empty", "zeros", "eye") \
                    and isinstance(st.value.func.value, ast.Name) \
                    and st.value.func.value.id == "np":
                shape, kind = self.m.allocs[tgt.id]
                cname = self.cn(tgt.id)
                if tgt.id not in self.outs_map:
                    dims = "".join(f"[{d}]" for d in shape)
                    self.lines.append(f"{indent}double {cname}{dims};")
                self.env[tgt.id] = shape
                if kind in ("zeros", "eye"):
                    total = 1
                    for d in shape:
                        total *= d
                    self.lines.append(
                        f"{indent}for (int _i = 0; _i < {total}; _i++)"
                        f" ((double *){cname})[_i] = 0.0;")
                if kind == "eye":
                    self.lines.append(
                        f"{indent}for (int _i = 0; _i < {shape[0]}; _i++)"
                        f" {cname}[_i][_i] = 1.0;")
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
                    if tgt.id in self.declared or tgt.id in self.env:
                        raise EmitError(f"{fn}: array local reassigned")
                    args = [self.expr(a) for a in st.value.args]
                    self.lines.append(
                        f"{indent}double {tgt.id}[{shape[0]}];")
                    self.lines.append(
                        f"{indent}{callee.c_name}"
                        f"({', '.join(args + [tgt.id])});")
                    self.env[tgt.id] = shape
                    return
                # array expression (matmul / transpose / elementwise)
                shp = _ashape(st.value, self.env, self.m.globals_)
                if shp is not None:
                    # The expression is lowered element-by-element into the
                    # destination buffer, so if that buffer is ALSO read on
                    # the RHS (e.g. t = a @ t, m = m + m.T) the in-flight
                    # writes corrupt later reads. The kernel-call path guards
                    # this as "array local reassigned"; the only safe alias is
                    # a bare `out = out`, which lower_array short-circuits.
                    if not (isinstance(st.value, ast.Name)
                            and st.value.id == tgt.id):
                        for sub in ast.walk(st.value):
                            if isinstance(sub, ast.Name) and sub.id == tgt.id:
                                raise EmitError(
                                    f"{fn}: array expression assigns to "
                                    f"{tgt.id!r} while reading it "
                                    f"(self-aliasing miscompiles)")
                    if tgt.id in self.outs_map or tgt.id in self.env:
                        dest = self.cn(tgt.id)
                        if tgt.id in self.outs_map \
                                and tgt.id not in self.env:
                            self.env[tgt.id] = shp
                        self.lower_array(st.value, dest, indent)
                    else:
                        dims = "".join(f"[{d}]" for d in shp)
                        self.lines.append(
                            f"{indent}double {tgt.id}{dims};")
                        self.env[tgt.id] = shp
                        self.lower_array(st.value, tgt.id, indent)
                    return
                rhs = self.expr(st.value)
                if tgt.id in self.declared:
                    self.lines.append(f"{indent}{tgt.id} = {rhs};")
                else:
                    self.declared.add(tgt.id)
                    self.lines.append(f"{indent}double {tgt.id} = {rhs};")
                return
            if isinstance(tgt, ast.Subscript) and isinstance(tgt.value, ast.Name):
                lhs = self.cn(tgt.value.id) + self._index_c(tgt.slice,
                                                            tgt.value.id)
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
                lhs = self.cn(st.target.value.id) \
                    + self._index_c(st.target.slice, st.target.value.id)
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
            self.loop_bounds[var] = bound
            self.lines.append(
                f"{indent}for (int {var} = 0; {var} < {bound}; {var}++) {{")
            for s in st.body:
                self.stmt(s, indent + "    ")
            self.lines.append(f"{indent}}}")
            self.loop_vars.discard(var)
            self.loop_bounds.pop(var, None)
            if st.orelse:
                raise EmitError(f"{fn}: for-else unsupported")
            return
        if isinstance(st, ast.If):
            # only these condition forms render with the outer parentheses C
            # requires; a bare name/subscript/call would emit `if flag {`.
            if not isinstance(st.test, (ast.Compare, ast.BinOp,
                                        ast.UnaryOp, ast.IfExp)):
                raise EmitError(
                    f"{fn}: unsupported if-condition "
                    f"{type(st.test).__name__} (only comparison/arithmetic "
                    f"conditions emit valid parenthesized C)")
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
        # hex float literals: ACSL parses decimal literals as exact
        # REALS, so `1e-05 <= n` is UNPROVABLE for the nearest double
        # (which sits below the real) — the boundary must be stated in
        # the double's own bits
        lo, hi = float(iv.lo).hex(), float(iv.hi).hex()
        if name in meta.param_arrays:
            n = meta.param_arrays[name]
            if n == 0:
                # size resolved only from loop-var subscripts: the ACSL bound
                # would be a vacuous `0 <= i < 0` and the signature a
                # non-ISO-C99 `const double v[0]` — neither is sound.
                raise EmitError(
                    f"{meta.py_name}: @contract on {name!r} whose array size "
                    f"resolves to 0 (sized by a loop bound, not a constant "
                    f"subscript); declare its length with @shapes")
            req.append(f"requires \\forall integer i; 0 <= i < {n} ==> "
                       f"{lo} <= {name}[i] <= {hi};")
        else:
            shp = meta.param_shapes.get(name)
            if shp is not None and len(shp) > 1:
                # a scalar range on a 2-D param would emit
                # `requires lo <= m <= hi;` for `const double m[.][.]`,
                # malformed ACSL that aborts Frama-C.
                raise EmitError(
                    f"{meta.py_name}: scalar-range @contract on multi-"
                    f"dimensional array parameter {name!r} (would emit "
                    f"malformed ACSL); contract its elements individually")
            req.append(f"requires {lo} <= {name} <= {hi};")
        spec.append(f"[in, {iv.to_annotation()}] {name};")
    acsl = "/*@\n  " + "\n  ".join(req) + "\n*/\n"
    ann = "/* [spec] { " + " ".join(spec) + " } */\n"
    return acsl + ann


def _signature(meta: _FuncMeta) -> str:
    parts = []
    for p in meta.param_order:
        shp = meta.param_shapes.get(p)
        if shp:
            dims = "".join(f"[{d}]" for d in shp)
            parts.append(f"const double {p}{dims}")
        else:
            parts.append(f"double {p}")
    if meta.outs:
        if len(meta.outs) == 1:
            dims = "".join(f"[{d}]" for d in meta.outs[0][1])
            parts.append(f"double out{dims}")
        else:
            for i, (_n, s) in enumerate(meta.outs):
                dims = "".join(f"[{d}]" for d in s)
                parts.append(f"double out{i}{dims}")
        ret = "void"
    else:
        ret = "double"
    return f"{ret} {meta.c_name}({', '.join(parts) or 'void'})"


def analyze_all(
    funcs: list[Callable],
) -> tuple[dict[str, _FuncMeta], dict[str, ast.FunctionDef]]:
    """Analyze all kernels and propagate array-ness through cross-kernel
    calls: a parameter never subscripted locally but passed where a
    callee expects an array IS an array of that size (two passes reach
    the fixpoint for the call depths the subset allows)."""
    metas: dict[str, _FuncMeta] = {}
    nodes: dict[str, ast.FunctionDef] = {}
    for f in funcs:
        meta, node = _analyze(f)
        metas[meta.py_name] = meta
        nodes[meta.py_name] = node
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
                        meta.param_shapes[arg.id] = (
                            meta.param_arrays[arg.id],)
    return metas, nodes


def emit_module(funcs: list[Callable], name: str = "podium_kernels",
                correctly_rounded: bool = False) -> str:
    """Emit a self-contained C99 translation unit for the given kernels.

    Order matters only for readability; prototypes are emitted first so
    cross-kernel calls resolve regardless of definition order.

    correctly_rounded: emit CORE-MATH cr_sin/cr_cos instead of libm
    sin/cos (and include "coremath.h"), making the transcendental
    kernels bit-exact against the correctly-rounded oracle rather than
    ULP-tolerant. Build with third_party/core-math/{sin,cos}.c.
    """
    metas, nodes = analyze_all(funcs)

    header = ["#include <math.h>"]
    if correctly_rounded:
        header.append('#include "coremath.h"')
    out = [f"/* {name} — generated by podium.emit.cemit; DO NOT EDIT.",
           "   Compile: gcc -std=c99 -O2 -ffp-contract=off (SSE2/binary64);",
           "   bit-exactness vs CPython depends on those flags. */",
           *header, ""]
    for meta in metas.values():
        out.append(_signature(meta) + ";")
    out.append("")
    for py_name, meta in metas.items():
        em = _Emitter(meta, metas, correctly_rounded=correctly_rounded)
        for st in nodes[py_name].body:
            em.stmt(st)
        out.append(_acsl(meta) + _signature(meta) + " {")
        out.extend(em.lines)
        out.append("}")
        out.append("")
    return "\n".join(out)
