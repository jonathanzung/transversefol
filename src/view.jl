
using Blink
#=
@eval AtomShell begin
    function init(; debug = false)
        electron() # Check path exists
        p, dp = port(), port()
        debug && inspector(dp)
        dbg = debug ? "--debug=$dp" : []
        proc = (debug ? run_rdr : run)(
            `$(electron()) --no-sandbox $dbg $mainjs port $p`; wait=false)
        conn = try_connect(ip"127.0.0.1", p)
        shell = Electron(proc, conn)
        initcbs(shell)
        return shell
    end
end
=#

function quickview(tup::NamedTuple; longitudes=true, obstructions=(tup.bt.ncusps < 3), contact_structures=false, h2=true, font_size=30, save_html=true, save_png=false,png_width::Int=1920, png_height::Int=1080, png_scale::Real=1, targets=[], fillings=Tuple{Int,Int}[])
    isosig = tup.isosig
    index = VeeringCensus.index(isosig)

    #include("batch/$(isosig).txt")
	Eupper = tup.Eupper
	Elower = tup.Elower


    #=
	if isosig == "eLMkbcddddedde_2100"
        dummy_candidate=random_cand(tup.bt, 1, DOWN)
		for pt in [(-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
			push!(Eupper, (pt, dummy_candidate))
			push!(Elower, (map(x->-x,pt), dummy_candidate))
		end
	end
    =#


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
    Econstr_upper = Envelope{Upper,Rational{Int},Cand{DiscreteHomeo{Tuple{Int,Int}}}}()
    Econstr_lower = Envelope{Lower,Rational{Int},Cand{DiscreteHomeo{Tuple{Int,Int}}}}()
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
        push!(longitudeDF, (ss=slopes(L), l=l, namedtuple_full(sss)..., text=string((normchi=normalizedchi(L),ss=ss,weights=l)), nchi = normalizedchi(L)))

        if is_fiber(l,tup.prep.top_bot_pairs)
            #@assert connected_components(l,fans)==1
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



    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in long_slopes)), init=CLIP),
                           maximum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in long_slopes)), init=-CLIP)] for i in 1:ncusps_unfilled]

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

    if length(Elower.A) > 0 && length(Eupper.A) > 0
        addtraces!(p, _plotjs(Elower, Eupper, name="Foliation slopes")...)
    end

    if haskey(tup, :Elowerbound)
        addtraces!(p, _plotjs(tup.Elowerbound, tup.Eupperbound, name="bound")...)
    end

    if obstructions
        if length(Econstr_lower.A) > 0
            addtraces!(p, _plotjs(Econstr_lower, Envelope{Upper}([([CLIP for i in 1:ncusps_unfilled], nothing)]), color=OBSTRUCTION_COLOUR, name="obstructions")...)
        end
        if length(Econstr_upper.A) > 0
            addtraces!(p, _plotjs(Envelope{Lower}([([-CLIP for i in 1:ncusps_unfilled],nothing)]), Econstr_upper, color=OBSTRUCTION_COLOUR, name="obstructions")...)
        end
    end


    function clip_df(df)
        return df
        return subset(df, :x => x->abs.(x).<=CLIP, :y => y->abs.(y).<=CLIP)
    end

    if nrow(longitudeDF) > 0
        #add_trace!(p, _plotjs(tup.Elong, color=LONGITUDE_COLOUR))
        if longitudes
            add_trace!(p, PlotlyJS.scatter(clip_df(longitudeDF); plotting_directives()..., marker=attr(line=attr(width=0), size=(ncusps_unfilled<=2 ? 25 : 10) ./ log.(4 .- longitudeDF[!,:nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers", name="Fibration slopes"))
        end
        if obstructions && nrow(constrDF) > 0
            add_trace!(p, PlotlyJS.scatter(clip_df(constrDF); plotting_directives()..., marker=attr(color=OBSTRUCTION_COLOUR, size=(ncusps_unfilled<=2 ? 5 : 3)), text=:text, mode="markers", name="obstructions"))
        end
    else
        println("no longitudes")
    end



    if !isempty(targets)
        target_args = ncusps == 1 ? (x=[t[1] for t in targets], y=zeros(length(targets))) :
                      ncusps == 2 ? (x=[t[1] for t in targets], y=[t[2] for t in targets]) :
                                    (x=[t[1] for t in targets], y=[t[2] for t in targets], z=[t[3] for t in targets], type="scatter3d")
        add_trace!(p, PlotlyJS.scatter(; target_args..., mode="markers",
            marker=attr(color="black", size=15, symbol="x"), name="filling slope"))
    end

    #=
    crevices = PEnvelope()
    for y in [(Vector{T}(x).+0.01, dummy_candidate) for x in crevices_general(Econstr_lower)]
        push!(crevices, y)
    end
    add_trace!(p, _plotjs(crevices, color=OBSTRUCTION_COLOUR))
    =#

    #add_trace!(p, _plotjs(Econstr_all, color=OBSTRUCTION_COLOUR))


    #add_trace!(p, _plotjs(Econstr_upper, color=OBSTRUCTION_COLOUR))
    #add_trace!(p, _plotjs(Econstr_lower, color=OBSTRUCTION_COLOUR))

    if false
        randE = random_trials(bt, nsubdivide=3, ntrials=1000000)
        randE2 = PEnvelope()
        for (x,c) in randE.A
            push!(randE2, (approximant_all_slopes(c::Candidate; time=10000), c))
        end

        #add_trace!(p, _plotjs(randE, color=TAUT_COLOUR))
        add_trace!(p, _plotjs(randE2, color=TAUT_COLOUR, name="random sample"))
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
                        marker=attr(color=LONGITUDE_COLOUR, size=10), name="H₂ boundary slopes"))
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
                    line=attr(color=LONGITUDE_COLOUR), name="H₂ boundary slopes"))
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
                    showscale=false, showlegend=true, name="H₂ boundary slopes"))
            end
        end
    end

    if save_html
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index).html"))
    end
    if save_png
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index).png"), width=png_width, height=png_height, scale=png_scale)
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


            for (s,cand) in Iterators.flatten([tup.Elong.A, tup.Eupper.A, tup.Elower.A])
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



function quickview_snappy(tup::NamedTuple; longitudes=true, h2=true, font_size=30, save_html=true, save_png=false, png_width::Int=1920, png_height::Int=1080, png_scale::Real=1, targets=[])
    isosig = tup.isosig
    index = VeeringCensus.index(isosig)
    bt = tup.bt
    ncusps = bt.ncusps
    transforms = degen_to_snappy_basis_change(bt)

    snappy_pairs = Envelopes.basis_change(tup.Elower, tup.Eupper, transforms)

    # Möbius transformation for coordinate i applied to a rational slope x
    mobius(i, x) = let a=transforms[i][1,1], b=transforms[i][1,2],
                       c=transforms[i][2,1], d=transforms[i][2,2]
        (c + d*x) / (a + b*x)
    end
    transform_slope(sss) = [mobius(i, sss[i]) for i in 1:ncusps]

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

    function plotting_directives()
        ncusps <= 2 ? (x=:x, y=:y, type="scatter") : (x=:x, y=:y, z=:z, type="scatter3d")
    end

    longitudeDF = DataFrame()
    long_slopes = []

    for l in tup.longitudes
        if !is_primitive(l) || connected_components(l, tup.prep.fans) != 1
            continue
        end
        L = Longitude(bt, l)
        ss = slopes(L)
        sss_s = transform_slope(map(slope_to_rat, ss))
        any(isinf, sss_s) && continue

        push!(long_slopes, sss_s)
        push!(longitudeDF, (ss=ss, l=l, namedtuple(sss_s)..., text=string((normchi=normalizedchi(L), ss=ss, weights=l)), nchi=normalizedchi(L)))
    end

    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/5, [Float64(x[i]) for x in long_slopes]), init=Float64(CLIP)),
                    maximum(filter(r -> abs(r) < CLIP/5, [Float64(x[i]) for x in long_slopes]), init=Float64(-CLIP))] for i in 1:ncusps]
    paddings = [0.5 * (y - x) for (x, y) in long_ranges]
    trimmed_ranges = [[r[1] - p, r[2] + p] for (r, p) in zip(long_ranges, paddings)]

    axes = if ncusps == 1
        (xaxis=attr(range=trimmed_ranges[1], minallowed=-CLIP, maxallowed=CLIP),
         yaxis=attr(showticklabels=false, range=[-1, 1]))
    elseif ncusps == 2
        (xaxis=attr(range=trimmed_ranges[1], title="cusp 1 (SnaPPy)", minallowed=-CLIP, maxallowed=CLIP),
         yaxis=attr(range=trimmed_ranges[2], title="cusp 2 (SnaPPy)", minallowed=-CLIP, maxallowed=CLIP))
    elseif ncusps == 3
        attr(scene=(xaxis=attr(range=trimmed_ranges[1]),
                    yaxis=attr(range=trimmed_ranges[2]),
                    zaxis=attr(range=trimmed_ranges[3])))
    end

    data = VeeringCensus.lookup_row(index)
    layout = Layout(title=attr(text="#$(index)  $(data[:isosig])  $(data[:names])  (SnaPPy coords)", font=attr(size=font_size));
                    axes..., font=attr(size=font_size),
                    hovermode=(ncusps < 3 ? "closest" : false))

    p = PlotlyJS.plot(layout)

    for (Elower_s, Eupper_s) in snappy_pairs
        addtraces!(p, _plotjs(Elower_s, Eupper_s, name="Foliation slopes")...)
    end

    if nrow(longitudeDF) > 0 && longitudes
        add_trace!(p, PlotlyJS.scatter(longitudeDF; plotting_directives()..., marker=attr(line=attr(width=0), size=(ncusps <= 2 ? 25 : 10) ./ log.(4 .- longitudeDF[!, :nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers", name="Fibration slopes"))
    end

    if !isempty(targets)
        target_args = ncusps == 1 ? (x=Float64[t[1] for t in targets], y=zeros(length(targets))) :
                      ncusps == 2 ? (x=Float64[t[1] for t in targets], y=Float64[t[2] for t in targets]) :
                                    (x=Float64[t[1] for t in targets], y=Float64[t[2] for t in targets], z=Float64[t[3] for t in targets], type="scatter3d")
        add_trace!(p, PlotlyJS.scatter(; target_args..., mode="markers",
            marker=attr(color="black", size=15, symbol="x"), name="filling slope"))
    end

    if h2
        h2_gens = compute_H2_rel_boundary(tup.prep.fans, tup.prep.tet_faces, tup.prep.face_coorientations)
        b1 = length(h2_gens)
        if b1 >= 1
            svs = [multislope_vec(bt, g) for g in h2_gens]

            function h2_snappy_slope(c, sp, sq)
                a, b, cc, d = transforms[c][1,1], transforms[c][1,2], transforms[c][2,1], transforms[c][2,2]
                denom = a*sp + b*sq
                denom == 0 && return NaN
                s = (cc*sp + d*sq) / denom
                abs(s) > CLIP ? NaN : s
            end

            if b1 == 1
                sv = svs[1]
                ss = [h2_snappy_slope(c, sv[c][1], sv[c][2]) for c in 1:ncusps]
                if all(!isnan, ss)
                    h2_args = ncusps == 1 ? (x=[ss[1]], y=[0.0]) :
                              ncusps == 2 ? (x=[ss[1]], y=[ss[2]]) :
                                            (x=[ss[1]], y=[ss[2]], z=[ss[3]], type="scatter3d")
                    add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="markers",
                        marker=attr(color=H2_COLOUR, size=10), name="H₂ boundary slopes"))
                end
            elseif b1 == 2
                θs = range(0.0, 2*Float64(π), length=4001)[1:end-1]
                sp(c, θ) = cos(θ)*svs[1][c][1] + sin(θ)*svs[2][c][1]
                sq(c, θ) = cos(θ)*svs[1][c][2] + sin(θ)*svs[2][c][2]
                xs = [h2_snappy_slope(1, sp(1,θ), sq(1,θ)) for θ in θs]
                h2_args = if ncusps == 1
                    (x=xs, y=zeros(length(θs)))
                elseif ncusps == 2
                    (x=xs, y=[h2_snappy_slope(2, sp(2,θ), sq(2,θ)) for θ in θs])
                elseif ncusps == 3
                    (x=xs, y=[h2_snappy_slope(2, sp(2,θ), sq(2,θ)) for θ in θs],
                     z=[h2_snappy_slope(3, sp(3,θ), sq(3,θ)) for θ in θs], type="scatter3d")
                end
                add_trace!(p, PlotlyJS.scatter(; h2_args..., mode="lines",
                    line=attr(color=H2_COLOUR), name="H₂ boundary slopes"))
            elseif b1 >= 3 && ncusps == 3
                φs = range(-Float64(π)/2, Float64(π)/2, length=400)
                θs = range(0.0, 2*Float64(π), length=800)[1:end-1]
                sp_s(c, φ, θ) = cos(φ)*cos(θ)*svs[1][c][1] + cos(φ)*sin(θ)*svs[2][c][1] + sin(φ)*svs[3][c][1]
                sq_s(c, φ, θ) = cos(φ)*cos(θ)*svs[1][c][2] + cos(φ)*sin(θ)*svs[2][c][2] + sin(φ)*svs[3][c][2]
                xs = [h2_snappy_slope(1, sp_s(1,φ,θ), sq_s(1,φ,θ)) for φ in φs, θ in θs]
                ys = [h2_snappy_slope(2, sp_s(2,φ,θ), sq_s(2,φ,θ)) for φ in φs, θ in θs]
                zs = [h2_snappy_slope(3, sp_s(3,φ,θ), sq_s(3,φ,θ)) for φ in φs, θ in θs]
                add_trace!(p, PlotlyJS.surface(x=xs, y=ys, z=zs,
                    colorscale=[[0, H2_COLOUR], [1, H2_COLOUR]],
                    showscale=false, showlegend=true, name="H₂ boundary slopes"))
            end
        end
    end

    if save_html
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index)_snappy.html"))
    end
    if save_png
        PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index)_snappy.png"), width=png_width, height=png_height, scale=png_scale)
    end
    flush(stdout)
    return p
end

function quickview_snappy(isosig::String; kwargs...)
    quickview_snappy(load(isosig); kwargs...)
end

function quickview_snappy(i::Int, ncusps::Int; kwargs...)
    quickview_snappy(load(VeeringCensus.lookup(i, ncusps)); kwargs...)
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
