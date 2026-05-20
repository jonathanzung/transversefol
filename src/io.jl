using Serialization

const BATCH_DIR = "/home/jonathan/batch" #joinpath(@__DIR__, "..", "batch")
const CLUSTER_BATCH_DIR = "/home/jonathan/engaging_sshfs/transversefol/batch"
const PYTHON_PATH = ENV["TRANSVERSEFOL_PYTHON"]
const BATCH_DIR = ENV["TRANSVERSEFOL_CACHE_DIR"]

function _load_prep(isosig::String)
    python = PYTHON_PATH
    script = joinpath(@__DIR__, "..", "pysrc", "prepare.py")
    json_str = read(`$python $script $isosig`, String)
    raw = JSON.parse(json_str)

    to_track(e) = (Int(e[1]), Int(e[2]))

    fans_py      = eval(Meta.parse(raw["fans"]))
    fc_py        = eval(Meta.parse(raw["face_coorientations"]))
    poles_py     = eval(Meta.parse(raw["poles"]))
    rungs_py     = eval(Meta.parse(raw["rungs"]))
    alledges_py  = eval(Meta.parse(raw["alledges"]))
    top_bot_py   = eval(Meta.parse(raw["top_bot_pairs"]))
    merid_py     = eval(Meta.parse(raw["meridian_dict"]))
    long_py      = eval(Meta.parse(raw["longitude_dict"]))
    degen_py     = eval(Meta.parse(raw["degeneracy"]))
    tet_faces_py = eval(Meta.parse(raw["tet_faces"]))

    return (
        fans          = [([to_track(e) for e in item[1]], [to_track(e) for e in item[2]]) for item in fans_py],
        face_coorientations = OffsetArrays.Origin(0)(Int[x for x in fc_py]),
        poles         = [[[to_track(e) for e in ladder] for ladder in cusp] for cusp in poles_py],
        rungs         = [[[to_track(e) for e in ladder] for ladder in cusp] for cusp in rungs_py],
        alledges      = [[to_track(e) for e in cusp] for cusp in alledges_py],
        top_bot_pairs = [to_track(x) for x in top_bot_py],
        meridian_dict  = Dict{Track,Int}(to_track(x) => Int(y) for (x, y) in merid_py),
        longitude_dict = Dict{Track,Int}(to_track(x) => Int(y) for (x, y) in long_py),
        degeneracy    = [(Int(x[1]), Int(x[2])) for x in degen_py],
        # tet_faces: for each tet, 4 (triangle_index, sign) pairs. sign is in Regina convention;
        # multiply by face_coorientations[tri_idx] to convert to veering co-orientation convention.
        tet_faces     = [[(Int(e[1]), Int(e[2])) for e in tet] for tet in tet_faces_py],
    )
end

function latest_save(filename)
    CUTOFF = Dates.datetime2unix(Dates.DateTime(2026,03,29,00,00) + Hour(4))
	locations=[]
    push!(locations, joinpath(BATCH_DIR, filename))
    push!(locations, joinpath(CLUSTER_BATCH_DIR, filename))

    locations = sort(filter(f->mtime(f) > CUTOFF, filter(isfile, locations)), by=mtime)
	if length(locations)==0
		println("not found")
        return nothing
	else
		@info "loading from $(locations[end])" Dates.unix2datetime(mtime(locations[end]))-Hour(4) #show in Eastern time zone
        return locations[end]
	end
end

function load(i::Int; kwargs...)
    load(VeeringCensus.lookup(i); kwargs...)
end
function load(i::Int, ncusps::Int; kwargs...)
    load(VeeringCensus.lookup(i,ncusps); kwargs...)
end

function loadstat(isosig::String; refresh=false)
    path = latest_save("$(isosig)_stat.json")
    if path==nothing || refresh
        return Dict{String,Any}("isosig"=>isosig)
    else
        return JSON.parsefile(path)
    end
end

function savestat(d::Dict)
    open(joinpath(BATCH_DIR, "$(d["isosig"])_stat.json"), "w") do io
        JSON.print(io, d)
    end
end

function load(isosig::String; refresh=false, refresh_longs=false, max_weight=100)
    prep = _load_prep(isosig)
    path = latest_save("$(isosig).jls")

	if path == nothing || refresh
        bt=BoundaryTriangulation(prep.fans,
                                 prep.face_coorientations,
                                 prep.alledges,
                                 prep.poles,
                                 prep.rungs,
                                 prep.meridian_dict,
                                 prep.longitude_dict)

        function filternan(A::Vector{T}) where {T}
            tmp = collect(filter(x->!(any(isinf.(x[1]))), A))
            return tmp
        end

		Elong, longitudes = compute_longitudes(bt, prep.fans, prep.top_bot_pairs, prep.tet_faces, prep.face_coorientations; max_weight=max_weight)
        Eupper = Envelope{Upper}([(x,set_roundmode(c, DOWN)) for (x,c) in filternan(Elong.A)])
        Elower = Envelope{Lower}([(x,set_roundmode(c, UP)) for (x,c) in filternan(Elong.A)])

		tup = (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong, longitudes=longitudes, isosig=isosig, prep = prep)
        save(tup)
	else
        tup = (deserialize(path)..., isosig=isosig)
        if refresh_longs
        	Elong, longitudes = compute_longitudes(tup.bt, tup.prep.fans, tup.prep.top_bot_pairs, tup.prep.tet_faces, tup.prep.face_coorientations; max_weight=max_weight)
            tup = (tup..., Elong=Elong, longitudes=longitudes)
            save(tup)
        end
	end
	return tup
end

function save(tup) #always save locally
	serialize(joinpath(BATCH_DIR, "$(tup.isosig).jls"), tup)
end


function computestat(isosig::String, f; name=String(Symbol(f)))
    d=loadstat(isosig)
    d[name] = f(load(isosig))
    savestat(d)
end

function computestat(indices::Union{AbstractVector{Int},AbstractRange{Int}}, f; name=String(Symbol(f)))
    @threads for ind in indices
        computestat(VeeringCensus.lookup(ind), f)
    end
end
function computestat(indices::Union{AbstractVector{Int},AbstractRange{Int}},ncusps::Int, f; name=String(Symbol(f)))
    @threads for ind in indices
        computestat(VeeringCensus.lookup(ind, ncusps), f)
    end
end

function reapstat(indices)
    df = DataFrame()
    for i in indices
        push!(df, loadstat(VeeringCensus.lookup(i)), cols=:union)
    end
    return df
end
function reapstat(indices, ncusps::Int)
    df = DataFrame()
    for i in indices
        push!(df, loadstat(VeeringCensus.lookup(i,ncusps)), cols=:union)
    end
    return df
end

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
    tup = load(isosig)

    @show length(tup.Elower.A), length(tup.Eupper.A)
    open("pointdump.ma","w") do f

        mathematica_print(f,[[x for (x,y) in tup.Elower.A],
        [x for (x,y) in tup.Eupper.A],
        find_s2_longitudes(tup)])
    end

    Econstr_lower, Econstr_upper = obstructions(tup; isosig=isosig)

    open("udump.ma", "w") do f
        mathematica_print(f,[crevices_general(Econstr_upper, CLIP), crevices_general(Econstr_lower, CLIP)])
    end
end

function dump_jankins_neumann()
    pts = Vector{Float64}[]
    for m in 1:100
        for a in 1:m-1
            if gcd(a,m)==1
                #push!(pts, Float64[(m-a)/m, a/m, 1+1/m])
                #push!(pts, Float64[a/m, 1-1/m, 1+a/m])
                #push!(pts, Float64[1-1/m, a/m, 1+a/m])
                push!(pts, Float64[a/m, (m-a)/m, 1/m])
                push!(pts, Float64[(m-a)/m, 1/m, a/m])
                push!(pts, Float64[1/m, (m-a)/m, a/m])
            end
        end
    end
    open("pointdump_jn.ma","w") do f
        mathematica_print(f,[[(0,0,0)], pts,
                             [(0,0,0)]])
    end
end

