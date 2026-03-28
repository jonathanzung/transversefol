function minchi(tup)
    chi=-Inf
	for l in tup.longitudes
        chi = max(chi, normalizedchi(Longitude(tup.bt,l)))
	end
    return chi
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

function compute_longitudes!(tup; max_weight=50)
	Elong, longitudes = compute_longitudes(tup.bt, tup.prep.fans, tup.prep.top_bot_pairs, tup.prep.tet_faces, tup.prep.face_coorientations; max_weight=max_weight)

    for (_, cand) in Elong.A
        push!(tup.Eupper, set_roundmode(cand, DOWN))
        push!(tup.Elower, set_roundmode(cand, UP))
    end
    save((tup..., Elong=Elong, longitudes=longitudes))
end


function compute_longitudes(bt, fans, top_bot_pairs, tet_faces, face_coorientations; max_weight=20)
    Elong = PEnvelope()
	long_dict = DefaultDict(()->[])
	longitudes = []

    longitudes_iter = find_longitudes_hom2(fans, top_bot_pairs, tet_faces, face_coorientations; max_weight=max_weight)
    for l in longitudes_iter
		if any([sum(abs.(x))==0 for x in slopes(Longitude(bt,l))])
			continue
		end
		ss = [y//x for (x,y) in slopes(Longitude(bt,l))]
		push!(long_dict[ss], l)
	end

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


function regimen(E::Envelope, target::Vector{T}; verbose=false, n_restarts=1, radius=0.2) where {T<:Real}
	#ncusps = length(E.A[1][1])

	println("phase 1")
	E = try_improve(E; nsubdivide=0, iters=30000, time=1000, target=target, radius=radius, beta=500, n_restarts=n_restarts)
	flush(stdout)
	if isinteractive() && verbose
		add_trace!(p, _plotjs(E))
	end
	println("phase 2")
	E = try_improve(E; nsubdivide=0, iters=300000, time=2000, target=target, radius=radius, beta=800, n_restarts=n_restarts)
	flush(stdout)
	if isinteractive() && verbose
		add_trace!(p, _plotjs(E))
	end
	println("phase 3")
	E = try_improve(E; nsubdivide=0, iters=1000000, time=2000, target=target, radius=radius, beta=1600, n_restarts=n_restarts)
	flush(stdout)
	if isinteractive() && verbose
		add_trace!(p, _plotjs(E))
	end
	return E
end

function runjob_closed(closed_isosig::String; kwargs...)
    parts = split(closed_isosig, '_')
    @assert length(parts)==3

    isosig = parts[1] * "_" * parts[2]

	tup = load(isosig, nlongs=0)

    snappy_target = eval(Meta.parse(parts[3]))
    local_target = [slope_to_rat(x*Slope(y)) for (x,y) in zip(snappy_to_degen_basis_change(tup.bt),snappy_target)]

    #runjob(isosig; merge((;kwargs...,), (target=local_target,))...)
	#tup = load(isosig, nlongs=0)

    @show local_target
    println(Envelopes.inclosure(tup.Eupper, local_target) && Envelopes.inclosure(tup.Elower, local_target))
end

function runjob(i::Int; kwargs...)
    runjob(VeeringCensus.lookup(i); kwargs...)
end

function runjob(i::Int, ncusps::Int; kwargs...)
    runjob(VeeringCensus.lookup(i,ncusps); kwargs...)
end

function boundjob(isosig::String; rt=0, thickness=24, ex=false, reg=false, nlongs=50, target=:crevices, fromscratch=false, doprune=false, preprune=false, refresh=false, verbose=false, bound=false, bound2=0, bound3=0, try_exclude=nothing, n_restarts=1, radius=0.2)

    	tup = load(isosig, refresh=refresh)

    if bound
        Eupperbound = Envelope{Upper}(eltype(Eupper.A)[])
        Elowerbound = Envelope{Lower}(eltype(Elower.A)[])


        g1 = all_cands(tup.bt, UP)
        @threads :greedy for c in g1
            push!(Eupperbound, (exact_slope(c),c))
		end

        g2 = all_cands(tup.bt, DOWN)
        @threads :greedy for c in g2
            push!(Elowerbound, (exact_slope(c),c))
		end
        @show [x[1] for x in Eupperbound.A]
        @show [x[1] for x in Elowerbound.A]


        save((tup..., Eupperbound=Eupperbound, Elowerbound=Elowerbound))
    end

    if try_exclude != nothing
        cand = basic_cand(tup.bt, DOWN)
        #ch=Channel{Tuple{Cand,Int}}(Inf)
        ch = []
        push!(ch, (cand, 1))

        i=0
        while !isempty(ch)
            if i%1000==0
                @show (i,length(ch))
            end
            i+=1
            (c,cusp) = pop!(ch)
            cdown = c
            cup = set_roundmode(c,UP)
            x = exact_slope(cdown)
            y = exact_slope(cup)
            @assert all(x .<= y)
            if i%1000000==0
                display(draw(c,1))
                display(draw(c,2))
                @show x, y
            end

            if !(all(x .<= try_exclude .<= y))
                continue
            elseif x==y
                println("found")
                break
            elseif x[cusp]==y[cusp]
                push!(ch, (c, mod1(cusp+1, tup.bt.ncusps)))
            else
                s = State(0//1, c.bt.rungs[cusp][1][1])
                while true
                    l = splittings(c, s)
                    if length(l) == 1
                        c=l[1]
                        s = trace_forwards(s,c)
                    else
                        for csplit in l
                            push!(ch, (csplit, mod1(cusp+1, tup.bt.ncusps)))
                        end
                        @assert false
                        break
                    end
                end
            end
        end
        @show (i,length(ch))
        println("excluded")
    end

    if bound3 > 0
        cand = basic_cand(tup.bt, DOWN)
        ch=Channel{Tuple{Cand,Int}}(Inf)
        Elowertmp = Envelope{Lower, Any, Any}()
        Euppertmp = Envelope{Upper, Any, Any}()

        push!(ch, (cand, 1))

        @threads :greedy for (c, cusp) in Iterators.take(ch,bound3)
            cdown = c
            cup = set_roundmode(c,UP)
            x = exact_slope(cdown)
            y = exact_slope(cup)
            @assert all(x .<= y)

            push!(Euppertmp, (x, cdown))
            push!(Elowertmp, (y, cup))

            if x==y
                continue
            elseif x[cusp]==y[cusp]
                push!(ch, (c, mod1(cusp+1, tup.bt.ncusps)))
            elseif !(Envelopes.inclosure(Euppertmp, y) && Envelopes.inclosure(Elowertmp, x))
                s = State(0//1, c.bt.rungs[cusp][1][1])
                while true
                    l = splittings(c, s)
                    if length(l) == 1
                        c=l[1]
                        s = trace_forwards(s,c)
                    else
                        for csplit in l
                            push!(ch, (csplit, mod1(cusp+1, tup.bt.ncusps)))
                        end
                        break
                    end
                end
            end
        end
        Eupperbound = Envelope{Upper}(copy(Euppertmp.A))
        Elowerbound = Envelope{Lower}(copy(Elowertmp.A))

        @threads :greedy for (c, cusp) in ch
            cdown = c
            cup = set_roundmode(c,UP)

            push!(Eupperbound, (exact_slope(cup), cup))
            push!(Elowerbound, (exact_slope(cdown), cdown))
            if isempty(ch)
                close(ch)
            end
        end
        @show [x[1] for x in Eupperbound.A]
        @show [x[1] for x in Elowerbound.A]


        save((tup..., Elower=Elowertmp, Eupper=Euppertmp, Eupperbound=Eupperbound, Elowerbound=Elowerbound))
    end

    if bound2 > 0
        cand = basic_cand(tup.bt, DOWN)

        current = [cand]

        Elowertmp = Envelope{Lower, Any, Any}()
        Euppertmp = Envelope{Upper, Any, Any}()

        for it in 0:bound2
            @show it
            for j in 1:tup.bt.ncusps
                next = []
                for c in current
                    _s = State(0//1, c.bt.rungs[j][1][1])
                    s=_s
                    verify_low(s,c,it-1)
                    for x in 1:it
                        s = trace_forwards(s,c)
                    end

                    #let's try to extend this state
                    for cnew in splittings(c,s)
                        verify_low(_s,cnew,it)
                        push!(next, cnew)
                    end
                end
                current = next



                g1 = current
                g2 = [set_roundmode(c,UP) for c in current]

                @threads :greedy for c in g1
                    push!(Euppertmp, (exact_slope(c),c))
                end

                @threads :greedy for c in g2
                    push!(Elowertmp, (exact_slope(c),c))
                end

                @show length(current)
                filter!(c-> begin
                                     x = exact_slope(c)
                                     y = exact_slope(set_roundmode(c,UP))
                                     @assert all(x .<= y)
                                     !(Envelopes.inclosure(Euppertmp, y) && Envelopes.inclosure(Elowertmp, x))
                                 end, current)
                @show length(current)
            end
        end

        Eupperbound = Envelope{Upper}(copy(Euppertmp.A))
        Elowerbound = Envelope{Lower}(copy(Elowertmp.A))

        g2 = current
        g1 = [set_roundmode(c,UP) for c in current]

        @threads :greedy for c in g1
            push!(Eupperbound, (exact_slope(c),c))
		end

        @threads :greedy for c in g2
            push!(Elowerbound, (exact_slope(c),c))
		end
        @show [x[1] for x in Eupperbound.A]
        @show [x[1] for x in Elowerbound.A]


        save((tup..., Elower=Elowertmp, Eupper=Euppertmp, Eupperbound=Eupperbound, Elowerbound=Elowerbound))
    end



end

function add_extreme_candidates!(isosig::String)
	if ex
		Elower2, Eupper2 = extreme_candidates(tup.bt)

		for x in Elower2.A
			push!(Elower, x)
		end
		for x in Eupper2.A
			push!(Eupper, x)
		end
	end
end

const TRY = (
            pool_size=1,
            n_walks=1, 
            iters=100000, 
            betastart=400, 
            betafinish=20000, 
            radius=0.1
        )
const TRYHARD = (
            pool_size=5,
            n_walks=10, 
            iters=200000, 
            betastart=400, 
            betafinish=20000, 
            radius=0.1
        )


function runjob(isosig::String; target=:crevices, refresh=false, showplots=true, candidates=[:longitudes], thickness=24, max_targets=100, optimization_args...)
	tup = load(isosig, refresh=refresh)

    @info "running $(tup.isosig)" target
    index = VeeringCensus.index(isosig)

	ncusps = tup.bt.ncusps
	Eupper = tup.Eupper
	Elower = tup.Elower

    upper_cands = []
    lower_cands = []

    if :longitudes in candidates
        for (_,c) in (tup.Elong.A)
            push!(upper_cands, set_roundmode(c,DOWN))
            push!(lower_cands, set_roundmode(c,UP))
        end
    end
    if :envelope in candidates
        for (_,c) in (tup.Eupper.A)
            push!(upper_cands, set_roundmode(c,DOWN))
        end
        for (_,c) in (tup.Elower.A)
            push!(lower_cands, set_roundmode(c,UP))
        end
    end       
    if :random in candidates
        for _ in 1:1000
            push!(upper_cands, random_cand(tup.bt,thickness,DOWN))
            push!(lower_cands, random_cand(tup.bt,thickness, UP))
        end
    end

    if target == :crevices
        upper_targets = crevices_general(Eupper,CLIP)
        lower_targets = crevices_general(Elower,CLIP)
    else
        upper_targets = [target]
        lower_targets = [target]
    end

    if length(upper_targets) > max_targets
        upper_targets = shuffle(upper_targets)[1:max_targets]
    end
    if length(lower_targets) > max_targets
        lower_targets = shuffle(lower_targets)[1:max_targets]
    end

    try_improve!(Eupper, upper_cands; targets=upper_targets, optimization_args...)
    save((tup..., Eupper=Eupper, Elower=Elower))
    try_improve!(Elower, lower_cands; targets=lower_targets, optimization_args...)
    save((tup..., Eupper=Eupper, Elower=Elower))

	p=quickview((tup..., Eupper=Eupper, Elower=Elower), targets=vcat(upper_targets, lower_targets))
	if isinteractive() && showplots
		display(p)
	end
    PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(isosig).html"))
    PlotlyJS.savefig(p, joinpath(BATCH_DIR, "$(index).html"))
    #save((tup..., Eupper=Eupper, Elower=Elower))
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

function obstructions(tup::NamedTuple)
    isosig = tup.isosig
	#include("batch/$(isosig).txt")

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

        if is_fiber(l,tup.prep.top_bot_pairs)
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

    return Econstr_lower, Econstr_upper
end


function review()
	include("batch/2cusp_manifest.txt")
	missing_isosigs = []
	interesting_isosigs = []
	for i in 1:120
		try
			tup = deserialize(joinpath(CLUSTER_BATCH_DIR, "$(isosigs[i]).jls"))
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
		f=joinpath(CLUSTER_BATCH_DIR, "$(isosig).jls")
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
    runjob(1,2, reg=true, rt=0, nlongs=2, refresh=true)
    Profile.Allocs.clear()
    Profile.init(n = 10^7, delay = 0.01)
    #Profile.Allocs.@profile sample_rate=0.0001 runjob(1, 2, reg=true, rt=0, nlongs=15, refresh=true)
    #PProf.Allocs.pprof()
    @profile runjob(1,2, reg=true, rt=0, nlongs=20, refresh=true)
end

function aggregate_bounds(X)
	include("batch/2cusp_manifest.txt")
    D=DataFrame()

	for i in 1:100
		isosig = isosigs[i]
        include("batch/$(isosig).txt")
		f=joinpath(CLUSTER_BATCH_DIR, "$(isosig).jls")
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

function dump_extremal(isosig)
    tup = load(isosig)

    @show length(tup.Elower.A), length(tup.Eupper.A)
    for (slopes,cand) in tup.Elower.A
        println()
        @show slopes
        show_trace(cand)
    end
end
