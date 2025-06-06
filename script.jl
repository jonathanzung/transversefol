using Serialization
using DataFrames
using StatProfilerHTML
using Revise
using Dates
using Profile,PProf


includet("search.jl")
includet("find_surface.jl")
includet("plotting.jl")
includet("envelopes.jl")

includet("envelopes.jl")

include("batch/2cusp_manifest.txt")

function mathematica_print(f::IO, l::Union{Array,Tuple})
	print(f,"{")
	for i in l[1:end-1]
		mathematica_print(f,i)
		print(f,",")
        #=
        if !(typeof(i) <: Real)
            println(f)
        end
        =#
	end
	mathematica_print(f,l[end])
	print(f,"}")
end

function mathematica_print(f::IO,l::Real)
	print(f,Float64(l))
end

function dump_points(isosig)
    tup = load(isosig)#deserialize("/home/jonathan/Dropbox/jonathan/transversefol/batch/$(isosig).jls")

    @show length(tup.Elower.A), length(tup.Eupper.A)
    open("pointdump.ma","w") do f

        mathematica_print(f,[[x for (x,y) in tup.Elower.A],
        [x for (x,y) in tup.Eupper.A],
        find_s2_longitudes(tup)])
    end

    Econstr_lower, Econstr_upper = obstructions(tup; isosig=isosig)

    open("udump.ma", "w") do f
        mathematica_print(f,[crevices_general(Econstr_upper), crevices_general(Econstr_lower)])
    end
end

function find_s2_longitudes(tup)
	D=Set()

    #=
	global bt=BoundaryTriangulation(fans, face_coorientations, firstrungs, alledges, rungs)	
	global ncusps = length(bt.firstrungs)
	global Elong=PEnvelope()
	global longitudeDF = DataFrame()
	global long_dict = DefaultDict(()->[])
	global longitudes = []
    =#

	for l in tup.longitudes
		if any([sum(abs.(x))==0 for x in slopes(Longitude(tup.bt,l))])
			continue
		end
		ss=[y//x for (x,y) in slopes(Longitude(tup.bt,l))]
		if normalizedchi(Longitude(tup.bt,l))==2 && all(abs.(ss) .< 100)
			push!(D,ss)
		end
	end
    return collect(D)
end

function compute_longitudes(bt; nlongs=100)
    Elong = PEnvelope()
	long_dict = DefaultDict(()->[])
	longitudes = []

    @show bt.rungs

    compute_homology(fans, top_bot_pairs)

    hom_classes=[map(x->x[1], y[1]) for y in bt.rungs[1:end-1]]
    @show hom_classes

	#ch = for l in find_longitudes_random(fans)
	#ch = find_longitudes_iterative(fans,1000) 
    ch = find_longitudes_hom(fans, top_bot_pairs)
    for l in ch #find_longitudes_hom(fans,top_bot_pairs)

        #=
        meridian = [[0,1],[1,0],[0,1]]

        @show (sum(l),slopes(Longitude(bt,l)),
               [det(hcat(x,y)) for (x,y) in zip(meridian, slopes(Longitude(bt,l)))]
              )
        #use this to find the class which intersects each meridian exactly once.
        =#
		if any([sum(abs.(x))==0 for x in slopes(Longitude(bt,l))])
            @show "rejected"
			continue
		end
		ss = [y//x for (x,y) in slopes(Longitude(bt,l))]

		if !haskey(long_dict, ss)
			@show length(long_dict)
			#flush(stdout)
		end
		push!(long_dict[ss], l)
		if length(long_dict) >= nlongs
			break
		end
	end
    close(ch)

	for (ss, ls) in long_dict

        #=
        if 1//4 in ss || -1//4 in ss
            @show ss
            @show ls
        end
        =#
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		push!(Elong, (ss, c))
        #push!(Elong, c)
		push!(longitudes, ls[i])
	end
	return Elong, longitudes
end

function latest_save(isosig)
	locations=[]
	push!(locations, "/home/jonathan/Dropbox/jonathan/transversefol/batch/$(isosig).jls")
	push!(locations, "/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls")

	locations = sort(filter(isfile, locations), by=mtime)
	if length(locations)==0
		println("not found")
        return nothing
	else
		println("loading from $(locations[end])")
        println(Dates.unix2datetime(mtime(locations[end]))-Hour(4)) #show in Eastern time zone
        return locations[end]
	end
end

function load(isosig; refresh=false, nlongs=100)
	println("setting up $(isosig), requesting $(nlongs) longitudes")
	flush(stdout)
	if !isfile("batch/$(isosig).txt")
		println("batch/$(isosig).txt not found, preparing it now")
		flush(stdout)
		run(`python3 prepare.py $(isosig)`)
	end

	include("batch/$(isosig).txt")

    path = latest_save(isosig)

	if path == nothing || refresh
		bt=BoundaryTriangulation(fans, face_coorientations, firstrungs, alledges, rungs)	
		ncusps = length(bt.firstrungs)

		Elong, longitudes = compute_longitudes(bt; nlongs=nlongs)
		Eupper = Envelope{Upper}(copy(Elong.A))
        Elower = Envelope{Lower}([(x,set_roundmode(c, UP)) for (x,c) in Elong.A])

        #Eupper= Envelope{Upper,Float64,Cand{DiscreteHomeo}}()
        #Elower= Envelope{Lower,Float64,Cand{DiscreteHomeo}}()


		tup = (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong, longitudes=longitudes)
		serialize("batch/$(isosig).jls", tup)
	else
		tup = deserialize(path)
		if length(tup.longitudes) < nlongs
			Elong, longitudes = compute_longitudes(tup.bt; nlongs=nlongs)
            for (x,c) in Elong.A
                push!(tup.Eupper, (x,c))
                push!(tup.Elower, (x,c))
            end
			tup = (bt=tup.bt, Eupper=tup.Eupper, Elower=tup.Elower, Elong=Elong, longitudes=longitudes)
			serialize("batch/$(isosig).jls", tup)
		end
	end
	println("done setup")
	flush(stdout)
	return tup
end

function save(isosig, tup) #always save locally
	serialize("batch/$(isosig).jls", tup)
end

function regimen(E::Envelope, target::Vector{T}; verbose=false) where {T<:Real}
	#ncusps = length(E.A[1][1])
	
	
	println("phase 1")
	E = try_improve(E; nsubdivide=0, iters=30000, time=1000, target=target, radius=0.001)
	flush(stdout)
	if isinteractive() && verbose
		add_trace!(p, _plotjs(E))
	end
	println("phase 2")
	E = try_improve(E; nsubdivide=0, iters=300000, time=2000, target=target, radius=0.001, beta=800)
	flush(stdout)
	if isinteractive() && verbose
		add_trace!(p, _plotjs(E))
	end
	println("phase 3")
	#E = try_improve(E; nsubdivide=1, iters=1000000, time=2000, target=target, radius=0.001, beta=1600)
	flush(stdout)
	if isinteractive() && verbose
		add_trace!(p, _plotjs(E))
	end
	return E
end

function runjob(i::Int; kwargs...)
	include("batch/2cusp_manifest.txt")
	isosig=isosigs[i]
	runjob(isosig; index=i, kwargs...)
end

function runjob(isosig::String; rt=0, ex=false, reg=false, nlongs=100, target=nothing, fromscratch=false, doprune=false, preprune=false, refresh=false, fix=false, verbose=false, index=0)
	tup = load(isosig, nlongs=nlongs, refresh=refresh)
    #=
	global p=quickview(tup; isosig=isosig)
	if isinteractive() && verbose
		display(p)
	end
    =#

	ncusps = length(tup.bt.firstrungs)


	Eupper = tup.Eupper
	Elower = tup.Elower

	if preprune
		prune!(Eupper)
		prune!(Elower)
	end

	if reg
		Eupper = regimen(Eupper, [CLIP for i in 1:ncusps]; verbose=verbose)
        serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
        Elower = regimen(Elower, [-CLIP for i in 1:ncusps]; verbose=verbose)
        serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
	end

    if target == :gaps
        Econstr_lower, Econstr_upper = obstructions(tup; isosig=isosig)
        Econstr_lower::Envelope{Lower}
        Econstr_upper::Envelope{Upper}

        upper_goals = collect(filter(x->!inclosure(Eupper, x), crevices_general(Econstr_lower)))
        lower_goals = collect(filter(x->!inclosure(Elower, x), crevices_general(Econstr_upper)))
        @show upper_goals
        @show lower_goals

        shuffle!(upper_goals)
        shuffle!(lower_goals)

        Euppertmp = if fromscratch
                        Envelope{Upper}(copy(tup.Elong.A))
                    else
                        Eupper
                    end

        Elowertmp = if fromscratch
                        Envelope{Lower}([(x,set_roundmode(c, UP)) for (x,c) in tup.Elong.A])
                    else
                        Elower
                    end

        @threads for target in upper_goals
            if !inclosure(Eupper, target)
                @show target
                Etmp = regimen(Envelope{Upper}(copy(Euppertmp.A)), target)
                for x in Etmp.A
                    push!(Eupper, x)
                end
                lock(Eupper.L) do
                    lock(Elower.L) do
                        serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
                    end
                end
            end
        end
        @threads for target in lower_goals
            if !inclosure(Elower, target)
                @show target
                Etmp = regimen(Envelope{Lower}(copy(Elowertmp.A)), target)
                for x in Etmp.A
                    push!(Elower, x)
                end
                lock(Eupper.L) do
                    lock(Elower.L) do
                        serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
                    end
                end
            end
        end
        serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
    elseif target != nothing
		Etmp = regimen(Envelope{Upper}(copy(
											if fromscratch
												tup.Elong.A
											else
                                                Eupper.A
											end
											)), target)
		for x in Etmp.A
			push!(Eupper, x)
		end
	end
    

	if ex
		Elower2, Eupper2 = extreme_candidates(tup.bt)

		for x in Elower2.A
			push!(Elower, x)
		end
		for x in Eupper2.A
			push!(Eupper, x)
		end
	end

	if rt>0
		randE = random_trials(tup.bt,ntrials=rt,thickness=24, roundmode=DOWN)
        #todo: multithread this
		@threads for (x,c) in randE.A
            push!(Eupper, (x,c))
		end

		randE = random_trials(tup.bt,ntrials=rt,thickness=16, roundmode=UP)
		@threads for (x,c) in randE.A
            #push!(Elower, (x,c))
        end

        serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
	end

	if doprune
		prune!(Eupper)
		prune!(Elower)
	end

	p=quickview((bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes); isosig=isosig, index=index)
	if isinteractive() && verbose
		display(p)
	end
    PlotlyJS.savefig(p, "batch/$(isosig).html")
    PlotlyJS.savefig(p, "batch/$(index).html")
	serialize("batch/$(isosig).jls", (bt=tup.bt, Eupper=Eupper, Elower=Elower, Elong=tup.Elong, longitudes=tup.longitudes))
    println("done job")
	flush(stdout)
end

function bench(c::Cand)
    @time for i in 1:1000
        exact_slope(c)
    end
    @time for i in 1:1000
        slope(c)
    end
end

function bench()
    tup = load(isosigs[1])
	c=random_cand(tup.bt,32,UP)
    bench(c)
end

function viewladderpole(i::Int)
    include("batch/2cusp_manifest.txt")
    run(`evince batch/$(isosigs[i]).pdf`)
end

function quickview(i::Int)
	#try
		include("batch/2cusp_manifest.txt")
		isosig=isosigs[i]
		quickview(isosig; index=i)
	#catch e
	#	@show e
	#end
end


function quickview(isosig::String; index=0)
	quickview(load(isosig); isosig=isosig, index=index)
end

#=
function constr_dict()
	include("batch/2cusp_manifest.txt")

    D=Dict()

    for i in 1:100
        errored=false
        local tup
        try
            tup = load_isosig(isosigs[i])
            include("batch/$(isosigs[i]).txt")
        catch error
            @show error
            errored=true
        end
        if errored
            continue
        end

        for l in tup.longitudes
            if is_fiber(l,top_bot_pairs)
                L=Longitude(tup.bt,l)
                for (k,(s,info)) in enumerate(constraints(L))
                    if all(!isnan(x) for x in s) && all(!isinf(x) for x in s)
                        if !haskey(D,info)
                            D[info]=[]
                        end
                        push!(D[info], i)
                    end
                end
            else
            end
        end

    end
    return D
end
=#

function obstructions(tup::NamedTuple; isosig="")
	include("batch/$(isosig).txt")

    bt=tup.bt

    Econstr = Envelope{Eq, Rational{Int}, Nothing}()
    Econstr_upper = Envelope{Upper,Rational{Int},Nothing}()
    Econstr_lower = Envelope{Lower,Rational{Int},Nothing}()
    longitudeDF = DataFrame()
    constrDF = DataFrame()
    long_slopes = []

    for l in tup.longitudes
        c=longitude_to_candidate(bt,l)
        L=Longitude(bt,l)
        ss=slopes(L)
        sss = map(slope_to_rat, ss)


        #b1=bound(tup.Eupper, sss[1])
        #b2=bound(tup.Elower, sss[1])
        push!(long_slopes, sss)

        if is_fiber(l,top_bot_pairs)
            #@assert connected_components(l,fans)==1
            for (s,info) in constraints(L)
                if all(!isnan(x) for x in s) && all(!isinf(x) for x in s) && info.npunc==1# && info.interior_prong >= 2
                    push!(Econstr, (s,nothing))
                    if info.dir[2] == -1
                        push!(Econstr_upper, (s,nothing))
                    else
                        @assert info.dir[2] == 1
                        push!(Econstr_lower, (s,nothing))
                    end

                end
            end
        else
            #println("nonfiber: $(sss)")
        end
    end

    #=
	Econstr=[PEnvelope() for i in 1:4]

	for l in tup.longitudes
        c=longitude_to_candidate(tup.bt,l)
        L=Longitude(tup.bt,l)
        ss=slopes(L)
        sss = [y//x for (x,y) in ss]

        if is_fiber(l, top_bot_pairs) && connected_components(l,fans)==1
            for (k,(s,info)) in enumerate(constraints(L))
                if all(!isnan(x) for x in s) && all(!isinf(x) for x in s) && info.npunc==1
                    push!(Econstr[k], (s,c))
                end
            end
        end
	end

	Econstr_upper = Envelope{Upper}()
	Econstr_lower = Envelope{Lower}()

	push!(Econstr_lower, Econstr[1])
	push!(Econstr_lower, Econstr[3])

	push!(Econstr_upper, Econstr[2])
	push!(Econstr_upper, Econstr[4])
    =#

    return Econstr_lower, Econstr_upper
end

function quickview(tup::NamedTuple; isosig="", index=0)
    include("batch/$(isosig).txt")
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


    @show length(tup.longitudes)


	bt = tup.bt
	ncusps = length(bt.firstrungs)

    function namedtuple(slopes)
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

    function plotting_directives()
        if ncusps ==1 || ncusps == 2
            return (x=:x, y=:y, type="scatter")
        elseif ncusps==3
            return (x=:x, y=:y, z=:z, type="scatter3d")
        end
    end


    Econstr = PEnvelope()
    Econstr_upper = Envelope{Upper,Float64,Cand{DiscreteHomeo}}()
    Econstr_lower = Envelope{Lower,Float64,Cand{DiscreteHomeo}}()
    longitudeDF = DataFrame()
    constrDF = DataFrame()
    long_slopes = []
    for l in tup.longitudes
        c=longitude_to_candidate(bt,l)
        L=Longitude(bt,l)
        ss=slopes(L)
        sss = map(slope_to_rat, ss)


        #b1=bound(tup.Eupper, sss[1])
        #b2=bound(tup.Elower, sss[1])
        push!(long_slopes, sss)
        push!(longitudeDF, (ss=slopes(L), l=l,  namedtuple(sss)..., text=string((normchi=normalizedchi(L),ss=ss,weights=l)), nchi = normalizedchi(L)))

        if is_fiber(l,top_bot_pairs)
            #@assert connected_components(l,fans)==1
            for (s,info) in constraints(L)
                if all(!isnan(x) for x in s) && all(!isinf(x) for x in s) && info.npunc==1# && info.interior_prong >= 2
                    push!(Econstr, (s,c))
                    push!(constrDF, (namedtuple(s)...,text=string(info)))
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



    long_ranges = [[minimum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in long_slopes)), init=CLIP),
                           maximum(filter(r -> abs(r) < CLIP/5, collect(x[i] for x in long_slopes)), init=-CLIP)] for i in 1:ncusps]

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
    
    axes = if ncusps == 1
        (
        xaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[1]
        ),
        yaxis=attr(
            showticklabels=false,
            range=[-1,1]
           ))
    elseif ncusps == 2
        (
        xaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[1],
            title="cusp 1 slope"
        ),
        yaxis=attr(
            showticklabels=false,
            range=trimmed_ranges[2],
            title="cusp 2 slope"
           ))
    elseif ncusps == 3
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

    layout = Layout(title="#$(index)   $(isosig)"; axes..., legend=attr(font=attr(
      size= 16)))

    p=PlotlyJS.plot(layout)

    #contact structures
    #addtraces!(p, _plotjs(Elower, Envelope{Upper,Float64,Cand{DiscreteHomeo}}([([CLIP for i in 1:ncusps], dummy_candidate)]), color=NEG_CONTACT_COLOUR, name="negative contact structures")...)

    #addtraces!(p, _plotjs(Envelope{Lower,Float64,Cand{DiscreteHomeo}}([([-CLIP for i in 1:ncusps],dummy_candidate)]), Eupper, color=POS_CONTACT_COLOUR, name="positive contact structures")...)

    addtraces!(p, _plotjs(Elower, Eupper, name="Z (foliated region)")...)

    if length(Econstr_lower.A) > 0
        addtraces!(p, _plotjs(Econstr_lower, Envelope{Upper}([([CLIP for i in 1:ncusps], nothing)]), color=OBSTRUCTION_COLOUR, name="obstructions")...)
    end
    if length(Econstr_upper.A) > 0
        addtraces!(p, _plotjs(Envelope{Lower}([([-CLIP for i in 1:ncusps],nothing)]), Econstr_upper, color=OBSTRUCTION_COLOUR, name="obstructions")...)
    end


    function clip_df(df)
        return df
        return subset(df, :x => x->abs.(x).<=CLIP, :y => y->abs.(y).<=CLIP)
    end


    #add_trace!(p, _plotjs(tup.Elong, color=LONGITUDE_COLOUR))
    add_trace!(p, PlotlyJS.scatter(clip_df(longitudeDF); plotting_directives()..., marker=attr(line=attr(width=0), size=(ncusps<=2 ? 25 : 10) ./ log.(4 .- longitudeDF[!,:nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers", name="fibrations"))
    add_trace!(p, PlotlyJS.scatter(clip_df(constrDF); plotting_directives()..., marker=attr(color=OBSTRUCTION_COLOUR, size=(ncusps<=2 ? 5 : 3)), text=:text, mode="markers", name="obstructions")) 



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

	PlotlyJS.savefig(p, "batch/$(index).html")
	#serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong))
	flush(stdout)
	p

	#interesting example isosigs[63]
end

function review()
	include("batch/2cusp_manifest.txt")
	missing_isosigs = []
	interesting_isosigs = []
	for i in 1:120
		try
			tup = deserialize("/home/jonathan/engaging_sshfs/transversefol/batch/$(isosigs[i]).jls")
			L=length(tup.Elower) + length(tup.Eupper)
			println("$i total weight $(L)")
			if L > 2
				push!(interesting_isosigs,i)
			end
		catch error
			@show error
			println("$i missing")
			push!(missing_isosigs,i)
		end
	end
	@show missing_isosigs
	@show interesting_isosigs

end

function verify(isosig::String)
	tup = load(isosig)
	@threads for (x,c) in tup.Elower.A
		@show x, approximant_all_slopes(c)
	end

	@threads for (x,c) in tup.Eupper.A
		@show x, approximant_all_slopes(c)
	end
end


function slope_to_rat(x::AbstractVector{T}) where {T<:Union{Int, Rational}}
    return x[2]//x[1]
end

function slope_to_rat(x::AbstractVector{T}) where {T<:Real}
    return x[2]/x[1]

end

function slope_to_twist(x::AbstractVector{T}) where {T<:Union{Int, Rational}}
    return x[1]//x[2]

end
function slope_to_twist(x::AbstractVector{T}) where {T<:Real}
    return x[1]/x[2]

end

function flagbad(range)
    include("batch/2cusp_manifest.txt")
    for n in range
        isosig=isosigs[n]
        include("batch/$(isosig).txt")
		f="/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls"
        if isfile(f)
            try
                tup = deserialize(f)
                if length(tup.Elong.A) != length(tup.Eupper.A)
                    for l in sort(tup.longitudes, by=sum)
                        L=Longitude(tup.bt,l)

                        ss=slopes(L)
                        npunctures = [gcd(a,b) for (a,b) in ss]
                        if is_fiber(l, top_bot_pairs) && connected_components(l,fans)==1

                            all_constraints = constraints(L)


                            q,p = ss[1]
                            s,r = ss[2]

                            sss = (p//q, r//s)

                            for i in 1:2
                                b1=rationalize(max(bound(tup.Eupper, sss, i), sss[i]),tol=1e-7)
                                b2=rationalize(min(bound(tup.Elower, sss, i), sss[i]),tol=1e-7)

                                upper, info = all_constraints[(i==1) ? 3 : 1]
                                lower, info = all_constraints[(i==1) ? 4 : 2]

                                if info.npunc==1 && info.interior_prong >= 2
                                    if b2 < lower[i] || b1 > upper[i]
                                        @show n,((T(b2),T(b1)),(T(lower[i]),T(upper[i])))
                                        @show n
                                        @goto here
                                    end
                                end   
                            end
                        end
                    end
                end
            catch e
                @show e
                @show n
            end
		end
        @label here

    end

end

function run_profile()
    runjob(1, reg=true, rt=0, nlongs=2, refresh=true)
    Profile.Allocs.clear()
    Profile.init(n = 10^7, delay = 0.01)
    Profile.Allocs.@profile sample_rate=0.0001 runjob(1, reg=true, rt=0, nlongs=100, refresh=true)
    PProf.Allocs.pprof()
end

function aggregate_bounds(X)
	include("batch/2cusp_manifest.txt")
    D=DataFrame()

	for i in 1:100
		isosig = isosigs[i]
        include("batch/$(isosig).txt")
		f="/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls"
        if isfile(f)
            try
                tup = deserialize(f)
                if length(tup.Elong.A) != length(tup.Eupper.A)
                    for l in sort(tup.longitudes, by=sum)
                        L=Longitude(tup.bt,l)

                        ss=slopes(L)
                        npunc = [gcd(a,b) for (a,b) in ss]
                        if is_fiber(l, top_bot_pairs) && connected_components(l,fans)==1

                            all_constraints = constraints(L)

                            q,p = ss[1]
                            s,r = ss[2]

                            b1=rationalize(max(bound(tup.Eupper, p//q), r//s),tol=1e-7)
                            b2=rationalize(min(bound(tup.Elower, p//q), r//s),tol=1e-7)

                            constr,info = all_constraints[1]
                            fol_upperbound = slope_to_twist(info.A*[1,b1])
                            fol_lowerbound = slope_to_twist(info.A*[1,b2])
                            if info.npunc==1 && info.chi>=X
                                push!(D, (info..., fol_upperbound=fol_upperbound, fol_lowerbound=fol_lowerbound, degen = slope_to_twist(info.degen_slope), chiinv=1/info.chi, text=string(i)))
                            end
                            
                        end
                    end
                end
            catch e
            end
		end
	end

	p=PlotlyJS.plot()

	add_trace!(p, PlotlyJS.scatter(D,x=:degen, y=:fol_lowerbound, marker=attr(line=attr(width=0), size=25 ./ log.(3 .- D[!,:chi]), color=LONGITUDE_COLOUR), text=:text, mode="markers"))

    for i in 1:-X
        for j in 1:i
            add_hline!(p, j/i)
        end
    end
    add_trace!(p, PlotlyJS.scatter(;x=0:1, y=0:1, mode="lines"))

    return p
end

#isosig = "siddhi2"
#isosig = "eLMkbcddddedde_2100"
#isosig = "gvLQQcdeffeffffaafa_201102" #L6a5
L6a5 = "gvLQQcdeffeffffaafa_201102"
#isosig = "gLLAQcdecfffhsermws_122201"
#isosig = "fLLQcbecdeepuwsua_20102"
#isosig = "fLLQcbeddeehhbghh_01110"
#isosig = "challenge2"
#isosig = "gLLPQbefefefhhhhhha_011102"
#isosig = "gLLPQcdfefefuoaaauo_022110"
#isosig = "gLLPQcdfefefuoaaauo_022110"
#setup(isosig)

#=
isosig = "eLMkbcddddedde_2100"
setup(isosig)

Profile.init(n=10^7, delay=0.01)
@profile setup(isosig)
using ProfileView
if isinteractive()
	ProfileView.view()
end
=#



#add_trace!(p, _plotjs(L6a2E, fill=true))
#add_trace!(p, _plotjs(accurate_E,fill=true))
#scatter_envelope!(p, accurate_E)
#scatter!(p, [-2,-1,-1/2,0],[1/2,1/3,1/6,0])
#scatter!(p,[-2,-1,-1/2],[1/3,1/6,0])


#=
if false
	c=subdivide(subdivide(longitude_to_candidate(bt,longitude)))
	println(all_slopes(c, time=50000))
	#=
	Profile.init(n=10^7, delay=0.01)
	vals1, candidate = @profile annealing(c->objective(all_slopes(c)), c, x->jiggle(x,0.005), 900, 1800, 2000000, verbose=true)
	println(all_slopes(candidate, time=50000))
	println(objective(all_slopes(candidate, time=10000)))
	=#
end
=#


#=
Profile.init(n=10^7, delay=0.01)
vals2, candidate = @profile annealing(c->objective2(all_slopes(c)), c, x->jiggle(x,0.03), 3000, 3000, 100000)
println(all_slopes(candidate, time=10000))
=#

#=
Profile.init(n=10^7, delay=0.01)
@profile random_trials(bt)
=#

#=
xs=0:0.01:1
p=plot(xs, [[f(x) for x in xs] for (J,f) in candidate.d if !J.inv])
display(p)
=#


#siddhi's example
#longitude  = 1/36
#meridian = 0
#46 triangles => 69 edges => H_1 has rank 24 => genus 12
#So the most we should expect is (36,1) - (2*12-1,0) = (13,1)
#1/13
#Looks like we're getting (1/3,1/3)
#Question is whether we can get (1/3+eps, 1/13-eps)
#Records: (1/13, 4/11) (1/3,1/3)



#bojun's example
#longitude = 1/48,   -1/4
#genus = 14
#prediction that we can get to 1/21, -1/4
#But we're getting 1/18
#



#L6a2
#Records
#(-2,1/2)   (-1, 1/3)    (-1/2, 1/6)    (-1/3, 1/9)    (-1/6, 1/18)
#
#
#pts = [(-2,1/2),   (-1, 1/3) ,   (-1/2, 1/6),    (-1/3, 1/9),    (-1/6, 1/18)]

#(0,0) is triangulated with 4 ideal triangles, has at least 2 punctures
# V - 3/2 *4 + 4 = V - 6 + 4 = V - 2. So V=2 => torus, V=4 => sphere. It's either a twice punctured torus or a 4 times punctured sphere 
#
#
#



# m125
# Records
# [[0.5, -2.0], [0.014492753623188406, -1.3714285714285714], [-1.375, 0.06666666666666667], [0.25, -1.5], [-1.3703703703703705, 0.014492753623188406], [0.0136986301369863, -1.368421052631579], [-2.0, 0.5], [-1.4, 0.1111111111111111], [0.1111111111111111, -1.4], [0.07142857142857142, -1.375], [0.0, 0.0], [-1.5, 0.25]]




#=
dummy_candidate=random_candidate(bt,0)
L6a2E=Envelope()
push!(L6a2E, (T[0,0], dummy_candidate))
for pt in [(-4,1/2+0.000001), (-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
	push!(L6a2E, (T[pt...], dummy_candidate))
end
=#


#problematic examples
#11 - endpoint
#38 - disaster
#35 - endpoint
#49 - endpoint
#71 - 1-prong surgery
#83 - not fibered
#91 - endpoint
#95 - disaster
#26 - not fibered
#160
#194 - 1-prong surgery
#214
#244 - 1-prong surgery

function dump_extremal(isosig)
    tup = load(isosig)#deserialize("/home/jonathan/Dropbox/jonathan/transversefol/batch/$(isosig).jls")

    @show length(tup.Elower.A), length(tup.Eupper.A)
    for (slopes,cand) in tup.Elower.A
        println()
        @show slopes
        show_trace(cand)
    end
end

function dump_candidate(c)

end


#bad_examples = [38, 95, 160, 214, 278, 338, 356, 370, 406, 448, 453, 470, 473, 485]
#
#
#
#
#



#For L6a5:
#Know from magic.py that in snappy coordinates
#
#degeneracy = (0,-1), (0,1), (-1,1)
#fiber = (4,1), (1,-1), (-1,0)
#meridian = (-1,0), (0,1), (-1,1)
#
#So in my coordinates, we know that the meridian is (0, infty, infty)
#Now we want to change coordinates, so that 0->infty, and 
