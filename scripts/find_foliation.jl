"""
find_foliation.jl — Find a taut foliation on a closed 3-manifold.

Usage:
    julia find_foliation.jl "s137(5,4)"
    julia find_foliation.jl "s137(5,4)" combinatorial

Searches for pseudo-Anosov flows on the given SnapPy manifold name using find_pA,
then for each resulting veering triangulation runs TransverseFol to find a taut
foliation transverse to the flow.  A foliation on the closed manifold exists iff
the Dehn filling slope lies inside the foliation envelope.
"""

import Pkg
Pkg.activate(@__DIR__; io=devnull)

using TransverseFol
using TransverseFol.Envelopes
using JSON

const PYTHON       = ENV["TRANSVERSEFOL_PYTHON"]
const FIND_PA_SCRIPT = joinpath(@__DIR__, "find_pA_flows.py")

# ---------------------------------------------------------------------------
# Step 1: call Python to enumerate pA flows for the given closed manifold
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Step 1: call Python to enumerate pA flows for the given closed manifold
# ---------------------------------------------------------------------------

"""
    find_pA_flows(manifold_name; method, count, max_drill, maxlength, isometry) -> Vector

Without `isometry`: returns `Vector{String}` of `<isosig>_[<slopes>]`.
With    `isometry`: returns `Vector{PAFlow}` with isometry data included.
"""
function find_pA_flows(manifold_name::String;
                       method       = "combinatorial",
                       count        = 10,
                       max_drill    = 3,
                       maxlength    = 3,
                       max_segments = 6,
                       max_tets     = 20,
                       isometry     = false)
    flags = isometry ? `--isometry` : ``
    cmd = `$(PYTHON) $(FIND_PA_SCRIPT) $(manifold_name)
           --method $(method) --count $(count) --max_drill $(max_drill)
           --maxlength $(maxlength) --max_segments $(max_segments) --max_tets $(max_tets)
           $(flags)`
    lines = readlines(cmd)
    isometry || return lines
    return parse_pA_flow_json.(lines)
end

# ---------------------------------------------------------------------------
# Step 2: convert a SnapPy filling slope to the Veering rational multislope
# ---------------------------------------------------------------------------

function snappy_filling_to_veering_rationals(tup, snappy_target)
    d2s = TransverseFol.degen_to_snappy_basis_change(tup.bt)
    n   = length(d2s)
    B   = inv(Envelopes.BasisChange([(0,0) for _ in 1:n], d2s, collect(1:n)))
    slopes = B(Tuple{Int,Int}[tuple(y...) for y in snappy_target])
    return [TransverseFol.slope_to_rat(s) for s in slopes]
end

# ---------------------------------------------------------------------------
# Step 3: check whether a multislope lies inside the achieved envelope
# ---------------------------------------------------------------------------

function is_achievable(tup, local_target)
    return (length(tup.Eupper.A) > 0 && length(tup.Elower.A) > 0
            && Envelopes.ininterior(tup.Eupper, local_target)
            && Envelopes.ininterior(tup.Elower, local_target))
end

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

"""
    find_foliation(manifold_name; method, pa_kwargs, runjob_kwargs...)

For each pA flow found by find_pA, runs the TransverseFol envelope search on
the underlying cusped veering triangulation, then checks whether the Dehn
filling slope is achievable.  Returns the first achievable closed isosig, or
`nothing`.
"""
function find_foliation(manifold_name::String;
                        method = "combinatorial",
                        pa_kwargs = NamedTuple(),
                        runjob_kwargs...)

    println("=== Searching for taut foliations on $(manifold_name) ===\n")

    closed_isosigs = unique(find_pA_flows(manifold_name; method=method, pa_kwargs...))
    println("Found $(length(closed_isosigs)) pA flow(s).")

    for isosig in closed_isosigs
        println(isosig)
    end
    
    if isempty(closed_isosigs)
        println("No pA flows found.")
        return nothing
    end

    # Group by base veering triangulation.
    # Each closed isosig has the form "<isosig>_[<slopes>]" where the isosig
    # itself contains exactly one underscore, e.g. "cPcbbbdxm_10_[(2,-3)]".
    base_to_closed = Dict{String, Vector{String}}()
    for s in closed_isosigs
        idx = findlast("_[", s)
        base = s[1:idx.start-1]
        push!(get!(base_to_closed, base, String[]), s)
    end

    for (base_isosig, closed_list) in sort(collect(base_to_closed), by=x->length(x[1]))
        println("\n--- Base triangulation: $(base_isosig) ---")

        tup = TransverseFol.load(base_isosig)

        for closed_isosig in closed_list
            idx = findlast("_[", closed_isosig)
            slopes_str    = closed_isosig[idx.start+1:end]
            snappy_target = eval(Meta.parse(slopes_str))
            local_target  = snappy_filling_to_veering_rationals(tup, snappy_target)

            tup = runjob(tup; runjob_kwargs..., target=local_target, showplots=false)

            achievable    = is_achievable(tup, local_target)

            println("  Filling $(slopes_str): Veering target = $(local_target),  achievable = $(achievable)")

            if achievable
                println("\n*** FOUND: taut foliation for $(manifold_name) via $(closed_isosig) ***")
                if isinteractive()
                    quickview(closed_isosig) |> display
                end
                return closed_isosig
            end
        end
    end

    println("\nNo taut foliation found among the discovered pA flows.")
    return nothing
end

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if !isempty(ARGS)
    manifold_name = ARGS[1]
    method        = length(ARGS) >= 2 ? ARGS[2] : "geodesic"
    find_foliation(manifold_name; method=method)
end
