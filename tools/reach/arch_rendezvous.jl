# ARCH spacecraft-rendezvous reachability regression.
#
# Consumes the hybrid-automaton JSON exported by
# `python -m podium.guidance.arch` (see export_model) and re-proves the
# benchmark safety properties with ReachabilityAnalysis.jl:
#   - line-of-sight cone and velocity octagon in the "attempt" mode,
#   - target-box avoidance in the "aborting" mode.
# Algorithm settings mirror the ARCH-COMP JuliaReach repeatability
# package (BOX delta=0.04, lazy clustering, box-template intersections).
# Exit code 0 iff every property is PROVEN — this is the CI gate.
#
# Usage: julia --project=tools/reach tools/reach/arch_rendezvous.jl model.json

using ReachabilityAnalysis, SparseArrays, JSON
using ReachabilityAnalysis: HalfSpace

function build_system(spec)
    n = length(spec["state"])
    modes = []
    for m in spec["modes"]
        A = Matrix{Float64}(hcat([Vector{Float64}(r) for r in m["A"]]...)')
        b = Vector{Float64}(m["b"])
        hss = [HalfSpace(Vector{Float64}(h["a"]), Float64(h["b"]))
               for h in m["invariant"]]
        inv = isempty(hss) ? Universe(n) : HPolyhedron(hss)
        push!(modes, @system(x' = A * x + b, x ∈ inv))
    end

    automaton = GraphAutomaton(length(modes))
    resetmaps = []
    for (i, tr) in enumerate(spec["transitions"])
        add_transition!(automaton, tr["from"], tr["to"], i)
        hss = [HalfSpace(Vector{Float64}(h["a"]), Float64(h["b"]))
               for h in tr["guard"]]
        push!(resetmaps, ConstrainedIdentityMap(n, HPolyhedron(hss)))
    end

    H = HybridSystem(automaton, [modes...], [resetmaps...],
                     [AutonomousSwitching()])
    X0 = Hyperrectangle(Vector{Float64}(spec["initial"]["center"]),
                        Vector{Float64}(spec["initial"]["radius"]))
    return IVP(H, [(spec["initial"]["mode"], X0)])
end

function properties(spec)
    n = length(spec["state"])
    p = spec["properties"]
    tan30 = p["tan30"]
    cx, cy = p["v_octagon_cx"], p["v_octagon_cy"]

    cone = HPolyhedron([
        HalfSpace(sparsevec([1], [-1.0], n), 100.0),
        HalfSpace(sparsevec([1, 2], [tan30, -1.0], n), 0.0),
        HalfSpace(sparsevec([1, 2], [tan30, 1.0], n), 0.0)])
    octagon = HPolyhedron([
        HalfSpace(sparsevec([3], [-1.0], n), cx),
        HalfSpace(sparsevec([3], [1.0], n), cx),
        HalfSpace(sparsevec([4], [-1.0], n), cx),
        HalfSpace(sparsevec([4], [1.0], n), cx),
        HalfSpace(sparsevec([3, 4], [1.0, 1.0], n), cy + cx),
        HalfSpace(sparsevec([3, 4], [1.0, -1.0], n), cy + cx),
        HalfSpace(sparsevec([3, 4], [-1.0, 1.0], n), cy + cx),
        HalfSpace(sparsevec([3, 4], [-1.0, -1.0], n), cy + cx)])
    target = BallInf(zeros(2), Float64(p["target_half_width"]))
    return cone, octagon, target, p["attempt_mode"],
           something(p["abort_mode"], -1)
end

function verify(sol, cone, octagon, target, attempt_mode, abort_mode)
    ok_los, ok_vel, ok_avoid = true, true, true
    for idx in findall(L -> L == attempt_mode, location.(sol))
        for R in sol[idx]
            ok_los &= set(R) ⊆ cone
            ok_vel &= set(R) ⊆ octagon
        end
    end
    if abort_mode > 0
        for idx in findall(L -> L == abort_mode, location.(sol))
            for R in sol[idx]
                box = overapproximate(Projection(R, [1, 2]), Hyperrectangle)
                ok_avoid &= isdisjoint(box, target)
            end
        end
    end
    return ok_los, ok_vel, ok_avoid
end

function main(path)
    spec = JSON.parsefile(path)
    println("model: ", spec["name"], "  abort_time=", spec["abort_time"])
    prob = build_system(spec)
    cone, octagon, target, attempt_mode, abort_mode = properties(spec)

    boxdirs = CustomDirections(collect(BoxDirections{Float64,Vector{Float64}}(5)))
    nclusters = spec["abort_time"] >= 0 ? 3 : 1
    stats = @timed solve(prob;
                alg=BOX(δ=0.04),
                clustering_method=LazyClustering(nclusters),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan=(0.0, Float64(spec["horizon"])))
    sol = stats.value
    ok_los, ok_vel, ok_avoid = verify(sol, cone, octagon, target,
                                      attempt_mode, abort_mode)
    println("reach time: ", round(stats.time; digits=2), " s")
    println("line_of_sight:    ", ok_los ? "PROVEN" : "VIOLATED/UNKNOWN")
    println("velocity_octagon: ", ok_vel ? "PROVEN" : "VIOLATED/UNKNOWN")
    println("abort_avoidance:  ", ok_avoid ? "PROVEN" : "VIOLATED/UNKNOWN")
    return ok_los && ok_vel && ok_avoid
end

all_ok = true
for path in ARGS
    global all_ok &= main(path)
end
exit(all_ok ? 0 : 1)
