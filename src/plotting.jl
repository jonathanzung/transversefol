using PlotlyJS
using WebIO
using JSON

const CLIP=25

# ---------------------------------------------------------------------------
# PAFlow: isometry data from find_pA_flows.py --isometry
# ---------------------------------------------------------------------------

# isosig:     base veering isosig (e.g. "eLMkbcddddedde_2100")
# fillings:   SnaPPy filling slopes for each cusp of M (the closed-manifold slopes)
# cusp_images[j]: cusp of MM that M's cusp j maps to (1-indexed), nothing if that cusp is filled
# cusp_maps[j]:   2×2 matrix (M SnaPPy → MM SnaPPy) for cusp j, nothing if filled
struct PAFlow
    isosig       :: String
    basis_change :: Envelopes.BasisChange   # M-SnaPPy -> MM-SnaPPy
end

function parse_pA_flow_json(line::String)
    d = JSON.parse(line)
    full_isosig = d["isosig"]
    idx = findlast("_[", full_isosig)
    isosig = full_isosig[1:idx.start-1]

    bc_slice = Tuple{Int,Int}[tuple(s...) for s in d["slice"]]
    perm     = Vector{Int}(d["perm"]) .+ 1  # 0-indexed -> 1-indexed
    bc_mats  = Matrix{Int}[Int[m[r][c] for r in 1:2, c in 1:2] for m in d["cusp_maps"]]

    #perm = collect(1:length(perm)) |> reverse
    #bc_mats  = Matrix{Int}[rand_sl2(true) for m in d["cusp_maps"]]
    #bc_mats = Matrix{Int}[[1 0; 0 1] for m in d["cusp_maps"]]

    return PAFlow(isosig, Envelopes.BasisChange(bc_slice, bc_mats, perm))
end

function rand_sl2(x)
    A = Matrix{Int}(I, 2, 2)
    for _ in 1:4
        n = rand(-3:3)
        if rand(Bool)
            A = [1 n; 0 1] * A
        else
            A = [1 0; n 1] * A
        end
    end
    if rand(Bool)
        A = [-1 0; 0 -1] * A
    end
    if x
        A = [-1 0; 0 1] * A
    end
    return A
end

function degen_to_snappy_basis_change_obj(bt)
    d2s = degen_to_snappy_basis_change(bt)
    #@info "d2s" d2s
    #d2s = Matrix{Int}[rand_sl2(true) for i in 1:bt.ncusps]
    #d2s = Matrix{Int}[[1 0; 0 1] for i in 1:bt.ncusps]
    n = length(d2s)
    return Envelopes.BasisChange([(0,0) for _ in 1:n], d2s, collect(1:n))
end

function inbounds(pt)
	return all(abs.(pt) .<= CLIP)
end

function clip_pt(pt)
    return [clamp(x, -CLIP, CLIP) for x in pt]
end

function staircase(E::Envelope{Upper})
	pts = sort!([x[1] for x in E.A], by=x->x[1])
	pts = unique!(clip_pt.(pts))
	return sort(vcat(pts,[(pts[i][1], pts[i+1][2]) for i in 1:length(pts)-1]), by=x->(x[1],-x[2]))
end

function staircase(E::Envelope{Lower})
	pts = sort!([x[1] for x in E.A], by=x->x[1])
	pts = unique!(clip_pt.(pts))
	return sort(vcat(pts,[(pts[i+1][1], pts[i][2]) for i in 1:length(pts)-1]), by=x->(x[1],-x[2]))
end

function plotjs(A::Vector{Envelope})
	PlotlyJS.plot([_plotjs(E) for E in A])
end

function plotjs(E::Envelope)
	plotjs(Envelope[E])
end
const TAUT_COLOUR = "#00cc96"
const MULTI_COLOURS = ["#00cc96", "#B6E880", "#FF97FF", "#FECB52", "#19D3F3", "#FF6692", "#FFA15A", "#72B7B2", "#AB63FA"]
const POS_CONTACT_COLOUR = "#eeee00"
const NEG_CONTACT_COLOUR = "#00eeee"
const OBSTRUCTION_COLOUR ="#ef553b"
const LONGITUDE_COLOUR = "#636EFA"
const H2_COLOUR = "#FF61FF"

function _plotjs(E::Envelope{S}; color=TAUT_COLOUR, name="") where {S<:Union{Upper,Lower}}
	pts = [x[1] for x in E.A]
	dim = length(pts[1])

	if dim==2
        marker = if color==nothing
            attr()
        else
            attr(color=color)
        end
		all_pts = staircase(E)
		PlotlyJS.scatter(x=Float64[x[1] for x in all_pts],y=Float64[x[2] for x in all_pts], mode="lines", line=attr(color=color), name=name, legendgroup=name)
	elseif dim==3
        marker = if color==nothing
            attr(size=3)
        else
            attr(color=color,size=3)
        end
		PlotlyJS.scatter(x=Float64[x[1] for x in pts],y=Float64[x[2] for x in pts], z=Float64[x[3] for x in pts], mode="markers", type="scatter3d", marker=marker)
	else
		@assert false
	end
end

function _plotjs(E::Envelope{S,T}; color=nothing) where {S,T}
	pts = filter(inbounds, [x[1] for x in E.A])


	if length(pts)==0
		return PlotlyJS.scatter(x=Float64[],y=Float64[], mode="markers")
	end
	dim = length(pts[1])

	if dim==2
        _pts = filter(x->abs(x[1]) <= CLIP && abs(x[2]) <= CLIP, pts)
        marker = if color==nothing
            attr()
        else
            attr(color=color)
        end
		PlotlyJS.scatter(x=Float64[x[1] for x in _pts],y=Float64[x[2] for x in _pts], mode="markers", marker=marker)
	elseif dim==3
        marker = if color==nothing
            attr(size=3)
        else
            attr(color=color,size=3)
        end

		PlotlyJS.scatter(x=Float64[x[1] for x in pts],y=Float64[x[2] for x in pts], z=Float64[x[3] for x in pts], mode="markers", type="scatter3d", marker=marker)
	else
		@assert false
	end
end

function _plotjs(E1::Envelope{Lower}, E2::Envelope{Upper}; color=TAUT_COLOUR, name="")
    #@info "envelopes" [x for (x,c) in E1.A] [x for (x,c) in E2.A]
    if isempty(E1.A) || isempty(E2.A)
        return []
    end
	dim = length(E1.A[1][1])
    if dim==1
        @assert length(E1.A)<=1
        @assert length(E2.A)<=1

        if length(E1.A)==0 || length(E2.A)==0 || E1.A[1][1] > E2.A[1][1]
            return []
        else
            #return [PlotlyJS.scatter(x=[E1.A[1][1], E2.A[1][1]],y=[0,0], mode="lines", name=name, legendgroup = name, line=attr(color=color))]
            return [PlotlyJS.scatter(x=Float64[clamp(E1.A[1][1][1],-CLIP,CLIP), clamp(E2.A[1][1][1],-CLIP,CLIP)],y=[0,0], mode="lines", name=name, legendgroup = name, line=attr(color=color, width=6))]
        end
    elseif dim==2
        rect_data = Envelopes.rectangles(E1, E2)
        rects = [PlotlyJS.scatter(
                    x=Float64[clamp(p1[1],-CLIP,CLIP), clamp(p2[1],-CLIP,CLIP), clamp(p2[1],-CLIP,CLIP), clamp(p1[1],-CLIP,CLIP), clamp(p1[1],-CLIP,CLIP)],
                    y=Float64[clamp(p1[2],-CLIP,CLIP), clamp(p1[2],-CLIP,CLIP), clamp(p2[2],-CLIP,CLIP), clamp(p2[2],-CLIP,CLIP), clamp(p1[2],-CLIP,CLIP)],
                    mode="none", fill="toself",
                    fillcolor=color * "33",
                    name=name, legendgroup=name, showlegend=(i==1))
                 for (i, (p1, p2)) in enumerate(rect_data)]

        all_pts1 = staircase(E1)
        all_pts2 = staircase(E2)

        staircases = [
            PlotlyJS.scatter(x=Float64[x[1] for x in all_pts1], y=Float64[x[2] for x in all_pts1], mode="lines", name=name, legendgroup=name, line=attr(color=color), showlegend=false),
            PlotlyJS.scatter(x=Float64[x[1] for x in all_pts2], y=Float64[x[2] for x in all_pts2], mode="lines", name=name, legendgroup=name, line=attr(color=color), showlegend=isempty(rects)),
        ]
        return vcat(rects, staircases)
	elseif dim==3
        cubes = [cube(p1, p2) for (p1, p2) in Envelopes.cuboids(E1, E2)]
        isempty(cubes) && return []
        return [combine_meshes(cubes; name=name)]
		#return [_plotjs(E1; color=color), _plotjs(E2, color=color)]
	end
end

function combine_meshes(v::Vector; name="") #a vector of tuples
    xs = collect(Iterators.flatten(m.x for m in v))
    ys = collect(Iterators.flatten(m.y for m in v))
    zs = collect(Iterators.flatten(m.z for m in v))
    facecolors = collect(Iterators.flatten(m.facecolor for m in v))

    is = Int[]
    js = Int[]
    ks = Int[]
    curroffset = 0

    for m in v
        append!(is, m.i .+ curroffset)
        append!(js, m.j .+ curroffset)
        append!(ks, m.k .+ curroffset)
        curroffset += length(m.x)
    end

    return mesh3d(x=xs,y=ys,z=zs,i=is,j=js,k=ks,facecolor=facecolors, flatshading=true, showlegend=true, name=name)
end

function cube(p1, p2)
    _p1 = Float64.(p1)
    _p2 = Float64.(p2)
    facecolor = repeat([
        "rgb(82,  188, 163)",  # teal       (front)
        "rgb(149, 105, 189)",  # violet     (bottom)
        "rgb(149, 105, 189)",  # violet     (top)
        "rgb(82,  188, 163)",  # teal       (back)
        "rgb(240, 128, 100)",  # coral      (left)
        "rgb(240, 128, 100)",  # coral      (right)
    ], inner=[2])
    t = (
               x=replace([0, 0, 1, 1, 0, 0, 1, 1], 0=> _p1[1], 1=>_p2[1]),
               y=replace([0, 1, 1, 0, 0, 1, 1, 0], 0=> _p1[2], 1=>_p2[2]),
               z=replace([0, 0, 0, 0, 1, 1, 1, 1], 0=> _p1[3], 1=>_p2[3]),
                i=[7, 0, 0, 0, 4, 4, 2, 6, 4, 0, 3, 7],
                j=[3, 4, 1, 2, 5, 6, 5, 5, 0, 1, 2, 2],
                k=[0, 7, 2, 3, 6, 7, 1, 2, 5, 5, 7, 6],
                facecolor=facecolor,
                flatshading=true

        )
    return t
end

function clear()
	deletetraces!(p,0:10)
end

function quickview(tup::NamedTuple; longitudes=true, obstructions=(tup.bt.ncusps < 3), contact_structures=false, h2=true, font_size=30, save_html=true, save_png=false,png_width::Int=1920, png_height::Int=1080, png_scale::Real=1, targets=[], fillings=Tuple{Int,Int}[])
    isosig = tup.isosig
    index = VeeringCensus.index(isosig)

	Eupper = tup.Eupper
	Elower = tup.Elower
	bt = tup.bt
    ncusps=bt.ncusps
    unfilled = isempty(fillings) ? collect(1:ncusps) : [i for i in 1:ncusps if fillings[i] == (0,0)]
    ncusps_unfilled = length(unfilled)
    if !isempty(fillings)
        Eupper = Envelopes.slice(Eupper, fillings)
        Elower = Envelopes.slice(Elower, fillings)
    end

    function namedtuple_full(slopes)
        @assert length(slopes)==ncusps
        if ncusps==1
            return (x=slopes[1], y=0)
        elseif ncusps==2
            return (x=slopes[1], y=slopes[2])
        elseif ncusps==3
            return (x=slopes[1], y=slopes[2], z=slopes[3])
        end
        @assert false
    end

    function namedtuple(slopes)
        @assert length(slopes)==ncusps_unfilled
        if ncusps_unfilled==1
            return (x=slopes[1], y=0)
        elseif ncusps_unfilled==2
            return (x=slopes[1], y=slopes[2])
        elseif ncusps_unfilled==3
            return (x=slopes[1], y=slopes[2], z=slopes[3])
        end
        @assert false
    end

    function plotting_directives()
        if ncusps_unfilled ==1 || ncusps_unfilled == 2
            return (x=:x, y=:y, type="scatter")
        elseif ncusps_unfilled==3
            return (x=:x, y=:y, z=:z, type="scatter3d")
        end
    end


    Econstr = PEnvelope()
    Econstr_upper = Envelope{Upper,Rational{Int},Nothing}()
    Econstr_lower = Envelope{Lower,Rational{Int},Nothing}()
    longitudeDF = DataFrame()
    constrDF = DataFrame()
    long_slopes = []
    for l in tup.longitudes
        if !is_primitive(l) || connected_components(l, tup.prep.fans) != 1
            continue
        end
        c=longitude_to_candidate(bt,l)
        L=Longitude(bt,l)
        ss=slopes(L)
        sss = map(slope_to_rat, ss)

        #b1=bound(tup.Eupper, sss[1])
        #b2=bound(tup.Elower, sss[1])
        push!(long_slopes, sss)
        push!(longitudeDF, (ss=slopes(L), l=l, namedtuple_full(sss)..., text=string((filledchi=normalizedchi(L),ss=Vector{Vector{Int}}(ss),weights=l)), nchi = normalizedchi(L)))

        if is_fiber(l,tup.prep.top_bot_pairs) || true
            #@assert connected_components(l,fans)==1
            
            #=
            for (s,info) in constraints(L)
                if all(!isnan(x) for x in s) && all(!isinf(x) for x in s) && info.npunc==1# && info.interior_prong >= 2
                    push!(Econstr, (s,c))
                    push!(constrDF, (namedtuple_full(s)..., s_raw=collect(s), text=string(info)))
                    if info.dir[2] == -1
                        push!(Econstr_upper, (s,c))
                    else
                        @assert info.dir[2] == 1
                        push!(Econstr_lower, (s,c))
                    end

                end
            end
            =#

            for s in constraints_conjecture_upper(L)
                push!(constrDF, namedtuple_full(s))
                push!(Econstr_lower, (s,nothing))
            end

            
            for s in constraints_conjecture_lower(L)
                push!(constrDF, namedtuple_full(s))
                push!(Econstr_upper, (s,nothing))
            end
            
        else
            #println("nonfiber: $(sss)")
        end
    end

    if !isempty(fillings)
        filling_compatible(sss) = all(fillings[i] == (0,0) || sss[i] == fillings[i][2]//fillings[i][1] for i in 1:ncusps)
        filling_compatible_float(sss) = all(fillings[i] == (0,0) || abs(sss[i] - fillings[i][2]/fillings[i][1]) < 1e-6 for i in 1:ncusps)

        longitudeDF_filtered = filter(row -> filling_compatible(map(slope_to_rat, row.ss)), longitudeDF)
        long_slopes = [[map(slope_to_rat, row.ss)[i] for i in unfilled] for row in eachrow(longitudeDF_filtered)]
        longitudeDF = DataFrame([merge(namedtuple([map(slope_to_rat, row.ss)[i] for i in unfilled]), (ss=row.ss, l=row.l, text=row.text, nchi=row.nchi)) for row in eachrow(longitudeDF_filtered)])

        constrDF_filtered = filter(row -> filling_compatible_float(row.s_raw), constrDF)
        constrDF = DataFrame([merge(namedtuple([row.s_raw[i] for i in unfilled]), (text=row.text,)) for row in eachrow(constrDF_filtered)])

        Econstr_upper = Envelopes.slice(Econstr_upper, fillings)
        Econstr_lower = Envelopes.slice(Econstr_lower, fillings)
    end



    envelope_slopes = vcat([Float64.(v) for (v,_) in Eupper.A], [Float64.(v) for (v,_) in Elower.A])
    all_range_slopes = vcat(long_slopes, envelope_slopes)
    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in all_range_slopes)), init=CLIP/5),
                           maximum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in all_range_slopes)), init=-CLIP/5)] for i in 1:ncusps_unfilled]

    paddings = [0.5 * (y-x) for (x,y) in long_ranges]
    trimmed_ranges = [[range[1]-pad, range[2]+pad] for (range,pad) in zip(long_ranges, paddings)]

    #=
    config = PlotConfig(modeBarButtonsToAdd=[
    "drawline",
    "drawopenpath",
    "drawclosedpath",
    "drawcircle",
    "drawrect",
    "eraseshape"
    ])
    =#
    config = PlotConfig(modeBarButtonsToAdd=[])

    #info = run(pipeline(`cat veering_census_with_data.txt`, `grep $(isosig)`))

    axes = if ncusps_unfilled == 1
        (
        xaxis=attr(
            showticklabels=true,
            range=trimmed_ranges[1],
            minallowed=-CLIP, maxallowed=CLIP
        ),
        yaxis=attr(
            showticklabels=false,
            range=[-1,1]
           ))
    elseif ncusps_unfilled == 2
        (
        xaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[1],
            title="cusp 1 slope",
            minallowed=-CLIP, maxallowed=CLIP
        ),
        yaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[2],
            title="cusp 2 slope",
            minallowed=-CLIP, maxallowed=CLIP
           ))
    elseif ncusps_unfilled == 3
        attr(scene=(
        xaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[1]
        ),
        yaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[2]
           ),

        zaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[3]
           )))

    end

    data = VeeringCensus.lookup_row(index)

    layout = Layout(title=attr(text="#$(index)   $(data[:isosig])   $(data[:depth])    $(data[:names])", font=attr(size=font_size));
                    axes...,
                    font=attr(size=font_size),
                    hovermode=(ncusps_unfilled < 3 ? "closest" : false))

    p=PlotlyJS.plot(layout)

    #contact structures
    if contact_structures
        dummy_candidate = random_cand(bt, 1, DOWN)
        addtraces!(p, _plotjs(Elower, Envelope{Upper,Rational{Int},Cand{DiscreteHomeo}}([(Rational{Int}[CLIP for i in 1:ncusps_unfilled], dummy_candidate)]), color=NEG_CONTACT_COLOUR, name="negative contact structures")...)

        addtraces!(p, _plotjs(Envelope{Lower,Rational{Int},Cand{DiscreteHomeo}}([(Rational{Int}[-CLIP for i in 1:ncusps_unfilled],dummy_candidate)]), Eupper, color=POS_CONTACT_COLOUR, name="positive contact structures")...)
    end

    if h2 && isempty(fillings)
        h2_gens = compute_H2_rel_boundary(tup.prep.fans, tup.prep.tet_faces, tup.prep.face_coorientations)
        b1 = length(h2_gens)
        if b1 >= 1
            svs = [multislope_vec(bt, g) for g in h2_gens]
            # Each svs[k] is a vector of ncusps 2-vectors [p, q]; slope ratio = q/p
            if b1 == 1
                sv = svs[1]
                # Show only if all denominators share the same sign (positive-denominator branch)
                if all(sv[c][1] > 0 for c in 1:ncusps) || all(sv[c][1] < 0 for c in 1:ncusps)
                    slope1(c) = sv[c][1] == 0 ? NaN : sv[c][2] / sv[c][1]
                    h2_args = ncusps == 1 ? (x=[slope1(1)], y=[0.0]) :
                              ncusps == 2 ? (x=[slope1(1)], y=[slope1(2)]) :
                                            (x=[slope1(1)], y=[slope1(2)], z=[slope1(3)], type="scatter3d")
                    add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="markers",
                        marker=attr(color=LONGITUDE_COLOUR, size=10), name="∂H₂ slopes"))
                end
            elseif b1 == 2
                # Iterate over [0, 2π); each projective point appears twice, but the
                # positive-denominator condition (all sp_c > 0) selects exactly one
                # representative per projective point, giving a clean connected arc.
                θs = range(0.0, 2*Float64(π), length=4001)[1:end-1]
                function h2_curve_slope(c, θ)
                    sp = cos(θ) * svs[1][c][1] + sin(θ) * svs[2][c][1]
                    sq = cos(θ) * svs[1][c][2] + sin(θ) * svs[2][c][2]
                    s = sp <= 0 ? NaN : sq/sp
                    abs(s) > CLIP ? NaN : s
                end
                mask = [all(cos(θ)*svs[1][c][1] + sin(θ)*svs[2][c][1] > 0 for c in 1:ncusps) for θ in θs]
                nanif(v, i) = mask[i] ? v : NaN
                xs = [nanif(h2_curve_slope(1, θ), i) for (i, θ) in enumerate(θs)]
                h2_args = if ncusps == 1
                    (x=xs, y=zeros(length(θs)))
                elseif ncusps == 2
                    (x=xs, y=[nanif(h2_curve_slope(2, θ), i) for (i, θ) in enumerate(θs)])
                elseif ncusps == 3
                    (x=xs, y=[nanif(h2_curve_slope(2, θ), i) for (i, θ) in enumerate(θs)],
                     z=[nanif(h2_curve_slope(3, θ), i) for (i, θ) in enumerate(θs)], type="scatter3d")
                end
                add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="lines",
                    line=attr(color=LONGITUDE_COLOUR), name="∂H₂ slopes"))
            elseif b1 >= 3 && ncusps == 3
                # Parametric surface over S²; positive-denominator branch selects
                # the unique representative (a₁:a₂:a₃) with all denom_c > 0,
                # covering each projective point in ℝP² exactly once.
                φs = range(-Float64(π)/2, Float64(π)/2, length=400)
                θs = range(0.0, 2*Float64(π), length=800)[1:end-1]
                function h2_surf_slope(c, φ, θ)
                    sp = cos(φ)*cos(θ)*svs[1][c][1] + cos(φ)*sin(θ)*svs[2][c][1] + sin(φ)*svs[3][c][1]
                    sq = cos(φ)*cos(θ)*svs[1][c][2] + cos(φ)*sin(θ)*svs[2][c][2] + sin(φ)*svs[3][c][2]
                    if !all(cos(φ)*cos(θ)*svs[1][cc][1] + cos(φ)*sin(θ)*svs[2][cc][1] + sin(φ)*svs[3][cc][1] > 0 for cc in 1:ncusps)
                        return NaN
                    end
                    s = sp == 0 ? Inf : sq/sp
                    abs(s) > CLIP ? NaN : s
                end
                xs = [h2_surf_slope(1, φ, θ) for φ in φs, θ in θs]
                ys = [h2_surf_slope(2, φ, θ) for φ in φs, θ in θs]
                zs = [h2_surf_slope(3, φ, θ) for φ in φs, θ in θs]
                add_trace!(p, PlotlyJS.surface(x=xs, y=ys, z=zs,
                    colorscale=[[0, LONGITUDE_COLOUR], [1, LONGITUDE_COLOUR]],
                    showscale=false, showlegend=true, name="∂H₂ slopes"))
            end
        end
    end



    if haskey(tup, :Elowerbound)
        addtraces!(p, _plotjs(tup.Elowerbound, tup.Eupperbound, name="bound")...)
    end

    if obstructions
        if length(Econstr_lower.A) > 0
            addtraces!(p, _plotjs(Econstr_lower, Envelope{Upper}([([1//0 for i in 1:ncusps_unfilled], nothing)]), color=OBSTRUCTION_COLOUR, name="obstructions")...)
        end
        if length(Econstr_upper.A) > 0
            addtraces!(p, _plotjs(Envelope{Lower}([([-1//0 for i in 1:ncusps_unfilled],nothing)]), Econstr_upper, color=OBSTRUCTION_COLOUR, name="obstructions")...)
        end
        cnstr_pts = vcat([(Float64.(v), :upper) for (v,_) in Econstr_upper.A],
                         [(Float64.(v), :lower) for (v,_) in Econstr_lower.A])
        if !isempty(cnstr_pts)
            cnstr_args = if ncusps_unfilled == 1
                (x=[v[1] for (v,_) in cnstr_pts], y=zeros(length(cnstr_pts)))
            elseif ncusps_unfilled == 2
                (x=[v[1] for (v,_) in cnstr_pts], y=[v[2] for (v,_) in cnstr_pts])
            else
                (x=[v[1] for (v,_) in cnstr_pts], y=[v[2] for (v,_) in cnstr_pts], z=[v[3] for (v,_) in cnstr_pts], type="scatter3d")
            end
            #=
            add_trace!(p, PlotlyJS.scatter(; cnstr_args..., mode="markers",
                marker=attr(color=OBSTRUCTION_COLOUR, size=(ncusps_unfilled<=2 ? 8 : 5), symbol=""),
                name="obstruction pts", legendgroup="obstructions", showlegend=false))
                =#
        end
    end

    if length(Elower.A) > 0 && length(Eupper.A) > 0
        addtraces!(p, _plotjs(Elower, Eupper, name="Foliation slopes")...)
    end


    function clip_df(df)
        return df
        return subset(df, :x => x->abs.(x).<=CLIP, :y => y->abs.(y).<=CLIP)
    end

    if nrow(longitudeDF) > 0
        #add_trace!(p, _plotjs(tup.Elong, color=LONGITUDE_COLOUR))
        if longitudes
            add_trace!(p, PlotlyJS.scatter(clip_df(longitudeDF); plotting_directives()..., marker=attr(line=attr(width=0), size=(ncusps_unfilled==1 ? 40 : ncusps_unfilled==2 ? 25 : 10) ./ log.(4 .- longitudeDF[!,:nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers", name="Fibration slopes"))
        end
        if obstructions && nrow(constrDF) > 0
            add_trace!(p, PlotlyJS.scatter(clip_df(constrDF); plotting_directives()..., marker=attr(color=OBSTRUCTION_COLOUR, size=(ncusps_unfilled==1 ? 10 : ncusps_unfilled==2 ? 5 : 3)), text=:text, mode="markers", legendgroup="obstructions",name="obstructions"))
        end
    else
        @info "no longitudes"
    end

    if !isempty(targets)
        target_args = ncusps == 1 ? (x=[t[1] for t in targets], y=zeros(length(targets))) :
                      ncusps == 2 ? (x=[t[1] for t in targets], y=[t[2] for t in targets]) :
                                    (x=[t[1] for t in targets], y=[t[2] for t in targets], z=[t[3] for t in targets], type="scatter3d")
        add_trace!(p, PlotlyJS.scatter(; target_args..., mode="markers",
            marker=attr(color="black", size=15, symbol="x"), name="targets"))
    end

    if save_html
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(tup.isosig).html"))
    end
    if save_png
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(tup.isosig).png"), width=png_width, height=png_height, scale=png_scale)
    end
	flush(stdout)

    on(p["click"]) do data
        @info data
        for point in data["points"]
            @info point
            if ncusps == 1
                coords = [point["x"]]
            elseif ncusps == 2
                coords = [point["x"],point["y"]]
            elseif ncusps == 3
                coords = [point["x"],point["y"],point["z"]]
            end

            if any(isnothing, coords)
                continue
            end

            @info (rationalize.(coords))

            for l in sort(tup.longitudes, by=sum)
                L=Longitude(tup.bt,l)

                ss=slopes(L)
                if rationalize.(coords)==map(slope_to_rat, ss)
                    @info (constraints(L))
                end
            end


            for (s,cand) in Iterators.flatten([tup.Eupper.A, tup.Elower.A, tup.Elong.A])
                if s == rationalize.(coords)
                    global lastcand = cand
                    slope_str = replace("($(join(rationalize.(coords), ",")))", "//" => "-")
                    for i in 1:ncusps
                        draw(cand, i; name="$(isosig)_$(slope_str)_cusp$(i)")
                    end
                    break
                end
            end
        end
    end

    return p

	#interesting example isosigs[63]
end


#=
function multiview(tups::Vector{<:NamedTuple}; longitudes=true, h2=true, font_size=30, save_html=true, save_png=false, png_width::Int=1920, png_height::Int=1080, png_scale::Real=1, targets=[], colors=MULTI_COLOURS)
    @assert !isempty(tups)
    ncusps = tups[1].bt.ncusps
    @assert all(t.bt.ncusps == ncusps for t in tups) "all tups must have the same ncusps"

    function plotting_directives()
        ncusps <= 2 ? (x=:x, y=:y, type="scatter") : (x=:x, y=:y, z=:z, type="scatter3d")
    end

    function namedtuple(slopes)
        if ncusps == 1
            return (x=Float64(slopes[1]), y=0.0)
        elseif ncusps == 2
            return (x=Float64(slopes[1]), y=Float64(slopes[2]))
        elseif ncusps == 3
            return (x=Float64(slopes[1]), y=Float64(slopes[2]), z=Float64(slopes[3]))
        end
        @assert false
    end

    # Collect all longitude slopes (from all tups) for axis range computation
    all_long_slopes = []
    for tup in tups
        bt = tup.bt
        B = degen_to_snappy_basis_change_obj(bt)
        for l in tup.longitudes
            !is_primitive(l) || connected_components(l, tup.prep.fans) != 1 && continue
            L = Longitude(bt, l)
            ss = slopes(L)
            result = B(Tuple{Int,Int}[(s[1], s[2]) for s in ss])
            result === nothing && continue
            any(r -> r[1] == 0, result) && continue
            push!(all_long_slopes, Float64[r[2]/r[1] for r in result])
        end
    end

    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/5, [Float64(x[i]) for x in all_long_slopes]), init=Float64(CLIP)),
                    maximum(filter(r -> abs(r) < CLIP/5, [Float64(x[i]) for x in all_long_slopes]), init=Float64(-CLIP))] for i in 1:ncusps]
    paddings = [0.5 * (y - x) for (x, y) in long_ranges]
    trimmed_ranges = [[r[1] - p, r[2] + p] for (r, p) in zip(long_ranges, paddings)]

    axes = if ncusps == 1
        (xaxis=attr(range=trimmed_ranges[1], minallowed=-CLIP, maxallowed=CLIP),
         yaxis=attr(showticklabels=false, range=[-1, 1]))
    elseif ncusps == 2
        (xaxis=attr(range=trimmed_ranges[1], title="cusp 1 (SnapPy)", minallowed=-CLIP, maxallowed=CLIP, showticklabels=false),
         yaxis=attr(range=trimmed_ranges[2], title="cusp 2 (SnapPy)", minallowed=-CLIP, maxallowed=CLIP, showticklabels=false))
    elseif ncusps == 3
        attr(scene=(xaxis=attr(range=trimmed_ranges[1]),
                    yaxis=attr(range=trimmed_ranges[2]),
                    zaxis=attr(range=trimmed_ranges[3])))
    end

    title_parts = ["#$(VeeringCensus.index(t.isosig)) $(t.isosig)" for t in tups]
    title_str = join(title_parts, "  |  ") * "  (SnapPy coords)"
    layout = Layout(title=attr(text=title_str, font=attr(size=font_size));
                    axes..., font=attr(size=font_size),
                    hovermode=(ncusps < 3 ? "closest" : false))

    p = PlotlyJS.plot(layout)

    for (k, tup) in enumerate(tups)
        bt = tup.bt
        index = VeeringCensus.index(tup.isosig)
        color = colors[mod1(k, length(colors))]
        B = degen_to_snappy_basis_change_obj(bt)

        label = length(tups) == 1 ? "Foliation slopes" : "Foliation slopes #$(index)"
        for (Elower_s, Eupper_s) in B * (tup.Elower, tup.Eupper)
            addtraces!(p, _plotjs(Elower_s, Eupper_s, color=color, name=label)...)
        end

        if longitudes
            longitudeDF = DataFrame()
            for l in tup.longitudes
                !is_primitive(l) || connected_components(l, tup.prep.fans) != 1 && continue
                L = Longitude(bt, l)
                ss = slopes(L)
                result = B(Tuple{Int,Int}[(s[1], s[2]) for s in ss])
                result === nothing && continue
                any(r -> r[1] == 0, result) && continue
                sss_s = Float64[r[2]/r[1] for r in result]
                push!(longitudeDF, (ss=ss, l=l, namedtuple(sss_s)..., text=string((normchi=normalizedchi(L), ss=ss, weights=l)), nchi=normalizedchi(L)))
            end
            if nrow(longitudeDF) > 0
                long_label = length(tups) == 1 ? "Fibration slopes" : "Fibration slopes #$(index)"
                add_trace!(p, PlotlyJS.scatter(longitudeDF; plotting_directives()..., marker=attr(line=attr(width=0), size=(ncusps <= 2 ? 25 : 10) ./ log.(4 .- longitudeDF[!, :nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers", name=long_label))
            end
        end

        if h2
            h2_gens = compute_H2_rel_boundary(tup.prep.fans, tup.prep.tet_faces, tup.prep.face_coorientations)
            b1 = length(h2_gens)
            if b1 >= 1
                svs_degen = [multislope_vec(bt, g) for g in h2_gens]
                # Pre-transform H2 generators into SnaPPy coordinates.
                # snappy_svs[k][c] = B.basis_change[c] * [p, q] in SnaPPy (p-coord, q-coord).
                snappy_svs = [[B.basis_change[c] * [svs_degen[k][c]...] for c in 1:ncusps] for k in 1:b1]

                snappy_slope(c, sp, sq) = sp == 0 ? NaN : (abs(sq/sp) > CLIP ? NaN : sq/sp)

                h2_label = length(tups) == 1 ? "H2 boundary slopes" : "H2 boundary slopes #$(index)"
                if b1 == 1
                    ss = [snappy_slope(c, snappy_svs[1][c][1], snappy_svs[1][c][2]) for c in 1:ncusps]
                    if all(!isnan, ss)
                        h2_args = ncusps == 1 ? (x=[ss[1]], y=[0.0]) :
                                  ncusps == 2 ? (x=[ss[1]], y=[ss[2]]) :
                                                (x=[ss[1]], y=[ss[2]], z=[ss[3]], type="scatter3d")
                        add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="markers",
                            marker=attr(color=LONGITUDE_COLOUR, size=10), name=h2_label))
                    end
                elseif b1 == 2
                    ts = range(0.0, 2*Float64(pi), length=4001)[1:end-1]
                    sp_i(c, t) = cos(t)*snappy_svs[1][c][1] + sin(t)*snappy_svs[2][c][1]
                    sq_i(c, t) = cos(t)*snappy_svs[1][c][2] + sin(t)*snappy_svs[2][c][2]
                    xs = [snappy_slope(1, sp_i(1,t), sq_i(1,t)) for t in ts]
                    h2_args = if ncusps == 1
                        (x=xs, y=zeros(length(ts)))
                    elseif ncusps == 2
                        (x=xs, y=[snappy_slope(2, sp_i(2,t), sq_i(2,t)) for t in ts])
                    elseif ncusps == 3
                        (x=xs, y=[snappy_slope(2, sp_i(2,t), sq_i(2,t)) for t in ts],
                         z=[snappy_slope(3, sp_i(3,t), sq_i(3,t)) for t in ts], type="scatter3d")
                    end
                    add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="lines",
                        line=attr(color=LONGITUDE_COLOUR), name=h2_label))
                elseif b1 >= 3 && ncusps == 3
                    phis = range(-Float64(pi)/2, Float64(pi)/2, length=400)
                    ts   = range(0.0, 2*Float64(pi), length=800)[1:end-1]
                    sp_i3(c, ph, t) = cos(ph)*cos(t)*snappy_svs[1][c][1] + cos(ph)*sin(t)*snappy_svs[2][c][1] + sin(ph)*snappy_svs[3][c][1]
                    sq_i3(c, ph, t) = cos(ph)*cos(t)*snappy_svs[1][c][2] + cos(ph)*sin(t)*snappy_svs[2][c][2] + sin(ph)*snappy_svs[3][c][2]
                    xs = [snappy_slope(1, sp_i3(1,ph,t), sq_i3(1,ph,t)) for ph in phis, t in ts]
                    ys = [snappy_slope(2, sp_i3(2,ph,t), sq_i3(2,ph,t)) for ph in phis, t in ts]
                    zs = [snappy_slope(3, sp_i3(3,ph,t), sq_i3(3,ph,t)) for ph in phis, t in ts]
                    add_trace!(p, PlotlyJS.surface(x=xs, y=ys, z=zs,
                        colorscale=[[0, LONGITUDE_COLOUR], [1, LONGITUDE_COLOUR]],
                        showscale=false, showlegend=true, name=h2_label))
                end
            end
        end
    end

    if !isempty(targets)
        target_args = ncusps == 1 ? (x=Float64[t[1] for t in targets], y=zeros(length(targets))) :
                      ncusps == 2 ? (x=Float64[t[1] for t in targets], y=Float64[t[2] for t in targets]) :
                                    (x=Float64[t[1] for t in targets], y=Float64[t[2] for t in targets], z=Float64[t[3] for t in targets], type="scatter3d")
        add_trace!(p, PlotlyJS.scatter(; target_args..., mode="markers",
            marker=attr(color="black", size=15, symbol="x"), name="filling slope"))
    end

    if save_html
        save_name = join([string(VeeringCensus.index(t.isosig)) for t in tups], "_") * "_snappy.html"
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, save_name))
    end
    if save_png
        save_name = join([string(VeeringCensus.index(t.isosig)) for t in tups], "_") * "_snappy.png"
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, save_name), width=png_width, height=png_height, scale=png_scale)
    end
    flush(stdout)
    return p
end
=#

function multiview(tup::NamedTuple; kwargs...)
    multiview([tup]; kwargs...)
end

function multiview(pa_flows::Vector{PAFlow}; kwargs...)
    multiview([(load(pf.isosig), pf) for pf in pa_flows]; kwargs...)
end

function multiview(pa_flow::PAFlow; kwargs...)
    multiview([pa_flow]; kwargs...)
end

"""
    multiview(isosigs_with_fillings; kwargs...)

Given a vector of strings of the form `"base_[(p,q),(0,0)]"`, calls
`isosig_isometries.py` to find isometries to a common manifold MM (determined
by the first entry), converts the results to PAFlow objects, loads the
corresponding veering triangulations, and plots everything in MM's SnaPPy
coordinate system.
"""
function multiview(isosigs_with_fillings::Vector{String}; MM=isosigs_with_fillings[1], kwargs...)
    @assert all(contains(s, "_[") for s in isosigs_with_fillings)

    python = "/home/jonathan/miniconda3/envs/sage/bin/python3"
    script = joinpath(@__DIR__, "..", "pysrc", "isosig_isometries.py")
    lines  = readlines(`$python $script $MM $isosigs_with_fillings`)
    pa_flows = parse_pA_flow_json.(lines)
    @info "flows" pa_flows
    multiview([(load(pf.isosig), pf) for pf in pa_flows]; kwargs...)
end

function multiview(isosig::String; kwargs...)
    multiview(load(isosig); kwargs...)
end

function multiview(i::Int, ncusps::Int; kwargs...)
    multiview(load(VeeringCensus.lookup(i, ncusps)); kwargs...)
end

"""
    multiview(tup_flows; title, longitudes, font_size, ...)

Plot foliation envelopes from multiple veering triangulations mapped into the
common SnaPPy coordinate system of the original manifold MM, using PAFlow
isometry data.  Each element of `tup_flows` is a `(tup, pa_flow)` pair.
"""
function multiview(tup_flows::Vector{<:Tuple{<:NamedTuple,PAFlow}};
                          title        = " ",
                          longitudes   = true,
                          obstructions = false,
                          h2           = true,
                          font_size    = 30,
                          save_html    = false,
                          save_png     = false,
                          png_width    :: Int   = 1920,
                          png_height   :: Int   = 1080,
                          png_scale    :: Real  = 1,
                          colors       = MULTI_COLOURS,
                          basis        :: Union{Envelopes.BasisChange, Nothing} = nothing,
                          flows        = eachindex(tup_flows),
                          LS_envelope :: Vector{<:Tuple{<:Envelope{Lower}, <:Envelope{Upper}}} = Tuple{Envelope{Lower}, Envelope{Upper}}[])
    tup_flows = tup_flows[collect(flows)]
    @assert !isempty(tup_flows)
    mm_ncusps = length(tup_flows[1][2].basis_change.perm)
    @assert all(length(pf.basis_change.perm) == mm_ncusps for (_, pf) in tup_flows)

    # Build reference inverse: maps MM-SnaPPy → M_1-degen (surviving cusps).
    # B_1 maps M_1-degen → MM-SnaPPy; strip its slice to get the bijective linear part.
    let (tup1, pf1) = tup_flows[1]
        B1 = pf1.basis_change * degen_to_snappy_basis_change_obj(tup1.bt)
        B1_bij = Envelopes.BasisChange(fill((0,0), mm_ncusps), B1.basis_change, B1.perm)
        global ref_inv = inv(B1_bij)
    end

    function display_B(tup, pa_flow)
        B_raw = pa_flow.basis_change * degen_to_snappy_basis_change_obj(tup.bt)
        B = ref_inv * B_raw
        return basis === nothing ? B : basis * B
    end

    function namedtuple(slopes)
        mm_ncusps == 1 ? (x=Float64(slopes[1]), y=0.0) :
        mm_ncusps == 2 ? (x=Float64(slopes[1]), y=Float64(slopes[2])) :
                         (x=Float64(slopes[1]), y=Float64(slopes[2]), z=Float64(slopes[3]))
    end
    plotting_directives() = mm_ncusps <= 2 ? (x=:x, y=:y, type="scatter") :
                                              (x=:x, y=:y, z=:z, type="scatter3d")

    # Apply BasisChange to a longitude slope vector; returns display Float64 slopes or nothing
    function apply_B_to_longitude(B, ss)
        result = B([(s[1], s[2]) for s in ss])
        result === nothing && return nothing
        disp_s = [Float64(r[2]) / Float64(r[1]) for r in result]
        any(isinf, disp_s) && return nothing
        return disp_s
    end

    # Collect all longitude slopes (for axis range computation)
    all_long_slopes = []
    for (tup, pa_flow) in tup_flows
        bt  = tup.bt
        B   = display_B(tup, pa_flow)
        for l in tup.longitudes
            !is_primitive(l) || connected_components(l, tup.prep.fans) != 1 && continue
            L  = Longitude(bt, l)
            ss = slopes(L)
            disp_s = apply_B_to_longitude(B, ss)
            disp_s === nothing && continue
            push!(all_long_slopes, disp_s)
        end
    end

    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/8, [x[i] for x in all_long_slopes]), init=Float64(CLIP)),
                    maximum(filter(r -> abs(r) < CLIP/8, [x[i] for x in all_long_slopes]), init=Float64(-CLIP))]
                   for i in 1:mm_ncusps]
    paddings = [0.25*(y-x) for (x,y) in long_ranges]
    trimmed_ranges = [[r[1]-p, r[2]+p] for (r,p) in zip(long_ranges, paddings)]

    axes = if mm_ncusps == 1
        (xaxis=attr(range=trimmed_ranges[1], minallowed=-CLIP, maxallowed=CLIP),
         yaxis=attr(showticklabels=false, range=[-1,1]))
    elseif mm_ncusps == 2
        (xaxis=attr(range=trimmed_ranges[1], title="cusp 1 slope", minallowed=-CLIP, maxallowed=CLIP, showticklabels=false),
         yaxis=attr(range=trimmed_ranges[2], title="cusp 2 slope", minallowed=-CLIP, maxallowed=CLIP, showticklabels=false))
    else
        attr(scene=(xaxis=attr(range=trimmed_ranges[1], showticklabels=false),
                    yaxis=attr(range=trimmed_ranges[2], showticklabels=false),
                    zaxis=attr(range=trimmed_ranges[3], showticklabels=false)))
    end

    function filled_isosig_str(pa_flow)
        fillings = pa_flow.basis_change.slice
        filled_slopes = [(p,q) for (p,q) in fillings]
        isempty(filled_slopes) && return pa_flow.isosig
        return pa_flow.isosig * "_" * string(filled_slopes)
    end

    title_str = isempty(title) ? "" :
        join(["$(k): $(filled_isosig_str(pf))" for (k, (_, pf)) in enumerate(tup_flows)], "   ")

    layout = Layout(title=attr(text=title_str, font=attr(size=font_size));
                    axes..., font=attr(size=font_size),
                    hovermode=(mm_ncusps < 3 ? "closest" : false))
    p = PlotlyJS.plot(layout)

    for (k, (tup, pa_flow)) in enumerate(tup_flows)
        bt  = tup.bt
        idx = k #VeeringCensus.index(pa_flow.isosig)
        color = colors[mod1(k, length(colors))]
        label = length(tup_flows) == 1 ? "Foliation slopes" : "Foliation slopes #$(idx)"

        B = display_B(tup, pa_flow)
        # Per-display-cusp matrix: output i uses surviving invperm[i]
        ip_B = invperm(B.perm)
        T_disp = [B.basis_change[ip_B[i]] for i in 1:mm_ncusps]

        for (El, Eu) in B * (tup.Elower, tup.Eupper)
            isempty(El.A) && isempty(Eu.A) && continue
            addtraces!(p, _plotjs(El, Eu, color=color, name=label)...)
        end

        if longitudes
            longitudeDF = DataFrame()
            for l in tup.longitudes
                !is_primitive(l) || connected_components(l, tup.prep.fans) != 1 && continue
                L  = Longitude(bt, l)
                ss = slopes(L)
                disp_s = apply_B_to_longitude(B, ss)
                disp_s === nothing && continue
                push!(longitudeDF, (namedtuple(disp_s)...,
                                    text=string((ss=slopes(L), weights=l)),
                                    nchi=normalizedchi(L)))
            end
            if nrow(longitudeDF) > 0
                long_label = length(tup_flows) == 1 ? "Fibration slopes" : "Fibration slopes #$(idx)"
                add_trace!(p, PlotlyJS.scatter(longitudeDF; plotting_directives()...,
                    marker=attr(line=attr(width=0),
                                size=(mm_ncusps <= 2 ? 25 : 10) ./ log.(4 .- longitudeDF[!,:nchi]),
                                color=LONGITUDE_COLOUR),
                    text=:text, mode="markers", name=long_label))
            end
        end

        if h2
            h2_gens = compute_H2_rel_boundary(tup.prep.fans, tup.prep.tet_faces, tup.prep.face_coorientations)
            b1 = length(h2_gens)
            if b1 >= 1
                svs = [multislope_vec(bt, g) for g in h2_gens]

                function h2_disp_slope(i, sp, sq)
                    a, b, c_coeff, d = T_disp[i][1,1], T_disp[i][1,2], T_disp[i][2,1], T_disp[i][2,2]
                    denom = a*sp + b*sq
                    denom == 0 && return NaN
                    s = (c_coeff*sp + d*sq) / denom
                    abs(s) > CLIP ? NaN : s
                end

                # unfilled_M[k] = M cusp index at position k among surviving coords of B
                unfilled_M = [j for j in 1:bt.ncusps if B.slice[j] == (0,0)]
                # M cusp for display coord i: surviving position is invperm[i]
                M_for_disp(i) = unfilled_M[ip_B[i]]

                h2_label = length(tup_flows) == 1 ? "∂H₂ slopes" : "∂H₂ slopes #$(idx)"
                if b1 == 1
                    sv_vals = [h2_disp_slope(i, svs[1][M_for_disp(i)][1], svs[1][M_for_disp(i)][2])
                               for i in 1:mm_ncusps]
                    if all(!isnan, sv_vals)
                        h2_args = mm_ncusps == 1 ? (x=[sv_vals[1]], y=[0.0]) :
                                  mm_ncusps == 2 ? (x=[sv_vals[1]], y=[sv_vals[2]]) :
                                                   (x=[sv_vals[1]], y=[sv_vals[2]], z=[sv_vals[3]], type="scatter3d")
                        add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="markers",
                            marker=attr(color=LONGITUDE_COLOUR, size=10), name=h2_label))
                    end
                elseif b1 == 2
                    ts = range(0.0, 2*Float64(pi), length=4001)[1:end-1]
                    sp_i(i, t) = cos(t)*svs[1][M_for_disp(i)][1] + sin(t)*svs[2][M_for_disp(i)][1]
                    sq_i(i, t) = cos(t)*svs[1][M_for_disp(i)][2] + sin(t)*svs[2][M_for_disp(i)][2]
                    xs = [h2_disp_slope(1, sp_i(1, t), sq_i(1, t)) for t in ts]
                    h2_args = if mm_ncusps == 1
                        (x=xs, y=zeros(length(ts)))
                    elseif mm_ncusps == 2
                        (x=xs, y=[h2_disp_slope(2, sp_i(2, t), sq_i(2, t)) for t in ts])
                    elseif mm_ncusps == 3
                        (x=xs, y=[h2_disp_slope(2, sp_i(2, t), sq_i(2, t)) for t in ts],
                         z=[h2_disp_slope(3, sp_i(3, t), sq_i(3, t)) for t in ts], type="scatter3d")
                    end
                    add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="lines",
                        line=attr(color=LONGITUDE_COLOUR), name=h2_label))
                elseif b1 >= 3 && mm_ncusps == 3
                    phis = range(-Float64(pi)/2, Float64(pi)/2, length=400)
                    ts   = range(0.0, 2*Float64(pi), length=800)[1:end-1]
                    sp_i3(i, ph, t) = cos(ph)*cos(t)*svs[1][M_for_disp(i)][1] + cos(ph)*sin(t)*svs[2][M_for_disp(i)][1] + sin(ph)*svs[3][M_for_disp(i)][1]
                    sq_i3(i, ph, t) = cos(ph)*cos(t)*svs[1][M_for_disp(i)][2] + cos(ph)*sin(t)*svs[2][M_for_disp(i)][2] + sin(ph)*svs[3][M_for_disp(i)][2]
                    xs = [h2_disp_slope(1, sp_i3(1, ph, t), sq_i3(1, ph, t)) for ph in phis, t in ts]
                    ys = [h2_disp_slope(2, sp_i3(2, ph, t), sq_i3(2, ph, t)) for ph in phis, t in ts]
                    zs = [h2_disp_slope(3, sp_i3(3, ph, t), sq_i3(3, ph, t)) for ph in phis, t in ts]
                    add_trace!(p, PlotlyJS.surface(x=xs, y=ys, z=zs,
                        colorscale=[[0, LONGITUDE_COLOUR], [1, LONGITUDE_COLOUR]],
                        showscale=false, showlegend=true, name=h2_label))
                end
            end
        end

        if obstructions
            # Compute obstruction envelopes in degen coords, then transform to display coords.
            Econstr_upper_degen = Envelope{Upper,Rational{Int},Nothing}()
            Econstr_lower_degen = Envelope{Lower,Rational{Int},Nothing}()
            for l in tup.longitudes
                !is_primitive(l) || connected_components(l, tup.prep.fans) != 1 && continue
                L = Longitude(bt, l)
                for s in constraints_conjecture_upper(L)
                    push!(Econstr_lower_degen, (s, nothing))
                end
                for s in constraints_conjecture_lower(L)
                    push!(Econstr_upper_degen, (s, nothing))
                end
            end
            obs_label = length(tup_flows) == 1 ? "obstructions" : "obstructions #$(idx)"
            ncusps_M = bt.ncusps
            Econstr_upper_inf = Envelope{Upper}([(Rational{Int}[1//0 for _ in 1:ncusps_M], nothing)])
            Econstr_lower_inf = Envelope{Lower}([(Rational{Int}[-1//0 for _ in 1:ncusps_M], nothing)])
            for (El_d, Eu_inf) in B * (Econstr_lower_degen, Econstr_upper_inf)
                isempty(El_d.A) || addtraces!(p, _plotjs(El_d, Eu_inf, color=OBSTRUCTION_COLOUR, name=obs_label)...)
            end
            for (El_inf, Eu_d) in B * (Econstr_lower_inf, Econstr_upper_degen)
                isempty(Eu_d.A) || addtraces!(p, _plotjs(El_inf, Eu_d, color=OBSTRUCTION_COLOUR, name=obs_label)...)
            end
            cnstr_pts = vcat(
                [(Float64.(v),) for (v,_) in (isempty(Econstr_upper_degen.A) ? [] : first(B * (Econstr_lower_inf, Econstr_upper_degen))[2].A)],
                [(Float64.(v),) for (v,_) in (isempty(Econstr_lower_degen.A) ? [] : first(B * (Econstr_lower_degen, Econstr_upper_inf))[1].A)]
            )
            if !isempty(cnstr_pts)
                cnstr_args = if mm_ncusps == 1
                    (x=[t[1][1] for t in cnstr_pts], y=zeros(length(cnstr_pts)))
                elseif mm_ncusps == 2
                    (x=[t[1][1] for t in cnstr_pts], y=[t[1][2] for t in cnstr_pts])
                else
                    (x=[t[1][1] for t in cnstr_pts], y=[t[1][2] for t in cnstr_pts],
                     z=[t[1][3] for t in cnstr_pts], type="scatter3d")
                end
                add_trace!(p, PlotlyJS.scatter(; cnstr_args..., mode="markers",
                    marker=attr(color=OBSTRUCTION_COLOUR, size=(mm_ncusps<=2 ? 8 : 5), symbol="x"),
                    name=obs_label, legendgroup=obs_label, showlegend=false))
            end
        end
    end

    for (El, Eu) in LS_envelope
        addtraces!(p, _plotjs(El, Eu, color=OBSTRUCTION_COLOUR,
                              name="L-space region",
                              )...)
    end

    if save_html
        save_name = join([string(VeeringCensus.index(pf.isosig)) for (_,pf) in tup_flows], "_") * "_MM_snappy.html"
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, save_name))
    end
    if save_png
        save_name = join([string(VeeringCensus.index(pf.isosig)) for (_,pf) in tup_flows], "_") * "_MM_snappy.png"
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, save_name), width=png_width, height=png_height, scale=png_scale)
    end
    flush(stdout)
    return p
end

function viewladderpole(i::Int)
    run(`evince $(joinpath(BATCH_DIR, "$(VeeringCensus.lookup(i)).pdf"))`)
end

function viewladderpole(i::Int, ncusps::Int)
    run(`evince $(joinpath(BATCH_DIR, "$(VeeringCensus.lookup(i,ncusps)).pdf"))`)
end

function viewladderpole(isosig::String)
    run(`evince $(joinpath(BATCH_DIR, "$(isosig).pdf"))`)
end

function quickview(i::Int, ncusps::Int; kwargs...)
    quickview(load(VeeringCensus.lookup(i, ncusps)); kwargs...)
end

function quickview(i::Int; kwargs...)
    quickview(load(VeeringCensus.lookup(i)); kwargs...)
end

function quickview(isosig::String; mode=:target, kwargs...)
    idx = findlast("_[", isosig)
    if isnothing(idx)
        quickview(load(isosig); kwargs...)
    else
        base_isosig  = isosig[1:idx.start-1]
        slopes_str   = isosig[idx.start+1:end]
        snappy_slopes = eval(Meta.parse(slopes_str))
        tup = load(base_isosig)
        basis = snappy_to_degen_basis_change(tup.bt)
        if mode == :slice
            fillings = [s == (0,0) ? (0,0) : begin v = basis[i] * Slope(s); (v[1], v[2]) end
                        for (i, s) in enumerate(snappy_slopes)]
            quickview(tup; fillings=fillings, kwargs...)
        else
            local_target = [slope_to_rat(basis[i] * Slope(snappy_slopes[i]))
                            for i in 1:length(snappy_slopes)]
            quickview(tup; targets=[local_target], kwargs...)
        end
    end
end
