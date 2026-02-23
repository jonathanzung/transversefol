using OffsetArrays
import Base: inv, getindex, setindex!, hash, push!, length, copy, show
using DataStructures
using Profile
#using ProfileView
using Base.Threads
using Measurements
using StaticArrays
using LinearAlgebra
using Envelopes
import Envelopes: Envelope
#import Plots: plot


const CLIP = 25 #don't worry about any slopes of absolute value bigger than CLIP



const Track = Tuple{Int,Int} #edge index and vertex number
const Slope = SVector{2,Int}
const MultiSlope{N} = SVector{N,Slope} where N #todo: make all vector of slopes into multislopes
abstract type Homeo end

#make this mutable, so that we can do simulated annealing by updating the same struct over and over again.

struct State{T}
	x::T#height in the edge
	e::Track#which track we're on
end

struct Junction
	index::Int
	inv::Bool
	left_len::Int
	right_len::Int
end

cartesianlength(::Type{Junction}) = 2 #cartesianlength tells the length of a type, as a tuple of ints, for the purposes of indexing.
cartesianlength(::Type{Tuple{S,T}}) where {S,T} = cartesianlength(S) + cartesianlength(T)
cartesianlength(::Type{Int}) = 1
tupleIndex(J::Junction) = (J.index+1, Int(J.inv)+1)
#tupleIndex((x,y)::Tuple{S,T}) where {S,T} = (tupleIndex(x)..., tupleIndex(y)...)
tupleIndex((x,y)::Tuple{Int,Int}) = (x+1,y+1)
tupleIndex(x::Int) = x+1
tupleIndex((x,J)::Tuple{Int,Junction}) = (x+1, J.index+1, Int(J.inv)+1)

struct ArrayDict{S,T,N}
	data::Array{T,N}
end

function copy(A::ArrayDict{S,T,N}) where {S,T,N}
    return ArrayDict{S,T,N}(copy(A.data))
end

function getindex(A::ArrayDict{S,T,N}, i::S) where {S,T,N}
	#@assert isassigned(A.data, tupleIndex(i)...)
	return A.data[CartesianIndex{N}(tupleIndex(i))]
end
function getindex(A::ArrayDict{S,T,2}, i::S) where {S,T}
	#@assert isassigned(A.data, tupleIndex(i)...)
	return A.data[tupleIndex(i)[1], tupleIndex(i)[2]]
end
function getindex(A::ArrayDict{S,T,3}, i::S) where {S,T}
	#@assert isassigned(A.data, tupleIndex(i)...)
	return A.data[tupleIndex(i)[1], tupleIndex(i)[2],tupleIndex(i)[3]]
end
function setindex!(A::ArrayDict{S,T,N}, val::T, i::S) where {S,T,N}
	A.data[CartesianIndex{N}(tupleIndex(i))]=val
end

function ArrayDict(D::Dict{S,T}) where {S,T}
	N=cartesianlength(S)
	indices = [tupleIndex(i) for i in keys(D)]
	ranges = (minimum(x[i] for x in indices):maximum(x[i] for x in indices) for i in 1:N)
	A=Array{T,N}(undef, 1 .+ length.(ranges)...)	
	#data = OffsetArray(A, ranges...)
	ret = ArrayDict{S,T,N}(A)

	for (k,v) in D
		ret[k]=v
	end
	return  ret
end

struct BoundaryTriangulation
	forward::ArrayDict{Track, Tuple{Int,Junction}, 2}
	backward::ArrayDict{Track, Tuple{Int,Junction}, 2}
	forwardfan::ArrayDict{Tuple{Int,Junction}, Track, 3} #leftward tracks
	backwardfan::ArrayDict{Tuple{Int,Junction}, Track, 3} #rightward tracks
	junctions::Vector{Junction}#doesn't include mappings for inverses

    poles::Vector{Vector{Vector{Track}}} #for each cusp, and for each ladder, the poles
	rungs::Vector{Vector{Vector{Track}}} #for each cusp, and for each ladder, the list of its rungs
	alledges::Vector{Vector{Track}}#edges appearing in each cusp

    zero_loops::Vector{Vector{Track}}#choice of slope zero in each cusp
	weights::ArrayDict{Track, Slope, 2}#these are stored relative to the taut coorientation, not regina's coorientation
	snappy_weights::ArrayDict{Track, Slope, 2}#these should also be stored relative to the taut coorientation. But these are dual weights, not weights

    ncusps::Int
    ntets::Int
	#fans is a list of pairs of fans
	
	function BoundaryTriangulation(fans, face_coorientations, alledges, poles, rungs, meridian_dict, longitude_dict)
		forward=Dict{Track, Tuple{Int,Junction}}()
		backward=Dict{Track, Tuple{Int,Junction}}()

		forwardfan=Dict{Tuple{Int,Junction}, Track}()
		backwardfan=Dict{Tuple{Int,Junction}, Track}()
		junctions=[]
		for (i,(f1,f2)) in enumerate(fans)
			J=Junction(i,false,length(f1),length(f2))
			push!(junctions,J)
			for (j,e) in enumerate(f1)
				_j=j-1
				coor = -face_coorientations[e[1]]
				forward[cwise(e,coor)]=(_j,J)
				forwardfan[(_j,J)] = cwise(e,coor)

				backward[ccwise(e,coor)]=(_j,inv(J))
				backwardfan[(_j,inv(J))] = ccwise(e,coor)
			end
			for (j,e) in enumerate(f2)
				coor = -face_coorientations[e[1]]
				_j=j-1
				forward[cwise(e,coor)]=(_j,inv(J))
				forwardfan[(_j,inv(J))] = cwise(e,coor)

				backward[ccwise(e,coor)]=(_j,J)
				backwardfan[(_j,J)] = ccwise(e,coor)
			end
		end
        ntets = length(junctions) #every tet has a canonical upper edge
        ncusps = length(rungs)

        snappy_weights = Dict{Track,Slope}()
        for l in alledges #iterate over cusps
            for track in l #for track in that cusp
                snappy_weights[track] = face_coorientations[track[1]] * Slope([meridian_dict[track], longitude_dict[track]])
            end
        end

		weights=Dict{Track,Slope}()
		for l in alledges
			for track in l
				weights[track] = Slope([0,0])
			end
		end


        zero_loops = Vector{Track}[]
        bt = new(ArrayDict(forward),
                 ArrayDict(backward),
                 ArrayDict(forwardfan),
                 ArrayDict(backwardfan),
                 junctions,

                 poles, 
                 rungs, 
                 alledges, 

                 zero_loops,
                 ArrayDict(weights),
                 ArrayDict(snappy_weights), 

                 ncusps, 
                 ntets)

		compute_weights!(bt)

		return bt
	end
end

struct Cand{S<:Homeo} #should it be mutable?
    bt::BoundaryTriangulation
    d::ArrayDict{Junction, S, 2}
end

function Base.show(io::IO, bt::BoundaryTriangulation)
    print(io, "BoundaryTriangulation(...,ncusps=$(bt.ncusps),ntets=$(bt.ntets))")
end
function Base.show(io::IO, c::Cand{S}) where {S}
    print(io, "Cand(...)")
end


includet("DiscreteHomeos.jl")

function Envelope()
    return Envelope{Upper,Float64,Cand{DiscreteHomeo{Tuple{Int,Int}}}}()
end

function PEnvelope()
    return Envelope{Eq,Float64,Cand{DiscreteHomeo{Tuple{Int,Int}}}}()
end

function right_tracks(bt::BoundaryTriangulation, J::Junction)
    return (bt.backwardfan[(i,J)] for i in 0:J.right_len-1)

end

function left_tracks(bt::BoundaryTriangulation, J::Junction)
    return (bt.forwardfan[(i,J)] for i in 0:J.left_len-1)
end

#=
function show(io::IO, c::Cand)
    println(io, "Cand")
    println(io, "exact slope: $(exact_slope(c))")
    println(io, "rung percentage: $(rung_percentage(c))")
end
=#

#=
struct Optimistic <: Homeo
	order::Vector{T,Bool} #(i,true) means the point at height i on the left
end

struct Pessimistic <: Homeo
end

function inv(o::S) where {S<: Union{Optimistic, Pessimistic}}
	return S([(x,!i) for (x,i) in o.order])
end
=#


#=
function (f::Optimistic)(y::T)
	for (h,b) in f.order
		if !b && 

end
=#

#=
function random_ordering(left_heights, right_heights)
    if length(left_heights)==2
        return Dir[BOTH, Dir[RIGHT for i in 1:length(right_heights)-2]..., BOTH]
    elseif length(right_heights) == 2
        return Dir[BOTH, Dir[LEFT for i in 1:length(left_heights)-2]..., BOTH]
    else
        k1 = rand(2:length(left_heights)-1)
        k2 = rand(2:length(right_heights)-1)

        return vcat(random_ordering(left_heights[1:k1], right_heights[1:k2])[1:end-1], random_ordering(left_heights[k1:end], right_heights[k2:end]))
    end
end
=#


function cwise((fnum,edge_index)::Tuple{Int,Int}, coor::Int)
	(fnum,mod(edge_index-coor,3))
end
function ccwise((fnum,edge_index)::Tuple{Int,Int}, coor::Int)
	(fnum,mod(edge_index+coor,3))
end


function Base.hash(J::Junction, h::UInt)
	return hash(J.index, hash(J.inv, h))
end



function inv(e::Junction)
	return Junction(e.index,!e.inv,e.right_len,e.left_len)
end

function getindex(c::Cand{H}, J::Junction) where {H}
	return c.d[J]
end

function setindex!(c::Cand{H}, f::H, J::Junction) where {H}
	c.d[J] = f
	c.d[inv(J)] = inv(f)
end

#=
function (c::Cand)(J::Junction, x::Rational{Int})
	@assert 0 <= x <= J.left_len
	return c[J](x // J.left_len) * J.right_len
end

function (c::Cand)(J::Junction, x)
	@assert 0 <= x <= J.left_len
	return c[J](x / J.left_len) * J.right_len
end
=#

function (c::Cand)(J::Junction, x)
    return c[J](x)
end

function set_roundmode(c::Cand{H}, r::RoundMode) where {H}
    cnew = Cand(c.bt, copy(c.d))
    for J in c.bt.junctions
        cnew[J] =set_roundmode(c[J], r)
        cnew[inv(J)] =inv(cnew[J])
    end
    return cnew
end

function subdivide(c::Cand{H}) where {H}
    cnew = Cand(c.bt, copy(c.d))
    for J in c.bt.junctions
        cnew[J] = subdivide(c[J])
        cnew[inv(J)] =inv(cnew[J])
    end
    return cnew
end

function set_roundmode(h::DiscreteHomeo, r::RoundMode)
    return DiscreteHomeo(h.ordering, h.dir, r)
end
function set_roundmode(h::DiscreteHomeo2, r::RoundMode)
    return DiscreteHomeo2(h.ordering_l, h.ordering_r, r)
end

#=
function subdivide(x::Vector{R}) where {R<: Rational}
    @assert length(x) >= 1
    y = R[]
    push!(y, x[1])
    for i in 2:length(x)
        push!(y, (y[end] + x[i])/2)
        push!(y,x[i])
    end
    return y
end
=#

function insert_left(h::DiscreteHomeo{T}, new) where {T}

end

function subdivide_left!(h::DiscreteHomeo{T}, old::T, new1::T, new2::T) where {T}
    r = searchsortedfirst(f.ordering, (old,old), by=a->a[f.dir])

    ind = rand(r)

    for i in r.start:ind
        h.ordering[i][f.dir] = new1
    end
    for i in ind:r.stop
        h.ordering[i][f.dir] = new2
    end
    insert!(f.ordering, ind, (new2, h.ordering[ind]))

    left = [x[1] for x in h.ordering]
    right = [x[2] for x in h.ordering]
    o=Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    for k in 2:length(h.ordering)
        push!(o, h.ordering[k])
        push!(o, h.ordering[k])
    end
    return DiscreteHomeo(subdivide(h.left_heights), subdivide(h.right_heights), o, h.dir, h.roundmode)
end

function copy(c::Cand)
    return Cand(c.bt, copy(c.d))
end

function jiggle(c::Cand{H}, r::T) where {H,T}
    cnew=Cand(c.bt,copy(c.d))
    for J in c.bt.junctions
        cnew[J] = jiggle(c[J], r)
        cnew[inv(J)] = inv(cnew[J])
    end
    return cnew
end

#=
function jiggle(f::DiscreteHomeo, r::T) where {T}
    dir = copy(f.ordering)


    for i in 1:3
        if rand() > 0.5
            valid_indices = filter(i-> dir[i]==BOTH || (dir[i] != dir[i+1] && dir[i+1] != BOTH), 2:length(dir)-1)
            index = rand(valid_indices)
            if dir[index]==BOTH
                if rand() > 0.5
                    dir = vcat(dir[1:index-1], [LEFT,RIGHT], dir[index+1:end])
                else
                    dir = vcat(dir[1:index-1], [RIGHT,LEFT], dir[index+1:end])
                end
            else
                if dir[index]==LEFT && dir[index+1]==RIGHT
                    dir = vcat(dir[1:index-1], Dir[BOTH], dir[index+2:end])
                elseif dir[index]==RIGHT && dir[index+1]==LEFT
                    dir = vcat(dir[1:index-1], Dir[BOTH], dir[index+2:end])
                end
            end
        end
    end
    return DiscreteHomeo(f.left_heights, f.right_heights, dir, f.dir, f.roundmode)
end
=#

#=
function jiggle(c::Candidate, r::T)
	cnew=Candidate(c.bt,Dict{Junction,Union{Piecewise,Linear}}())
	for J in c.bt.junctions
		#cnew[J] = (rand() < 0.3 ?  jiggle(c[J],r) : c[J])
		cnew[J] = jiggle(c[J],r)
	end
	return cnew
	#return a new candidate 
end

function jiggle(p::Piecewise, r::T)
	Piecewise(jiggle(p.left,r), jiggle(p.right,r), jiggle(p.x,r), jiggle(p.fx,r))
end

function jiggle(x::T, r::T)
	@assert 0 <= x <= 1
	#if rand() < 0.2
		return uniform(max(0, x-r), min(1,x+r))
	#else
	#	return x
	#end
end

function jiggle(x::Linear, r::T)
	Linear()
end
=#

function print_junction(bt::BoundaryTriangulation, J::Junction)
	println("Junction")
	for i in 0:J.left_len-1
		@show i, bt.forwardfan[(i,J)]
	end
	for i in 0:J.right_len-1
		@show i, bt.backwardfan[(i,J)]
	end
end

#junction_crossings=Dict()

#=
function trace_forwards(s::State{T}, c::Cand) where {T}
    #@show s.x
    #@show s.e
	@assert 0<= s.x <=1
	i,J = c.bt.forward[s.e]
	
#	print_junction(bt,J)
    fx = c(J,i+s.x)
    #@show fx


    j=Int(floor(fx))
	#@show j
	
	if j==J.right_len
        @assert fx==j
		#fx = j
		j=j-1
	end

    #println("next\n\n")

    return State{T}(fx-j, c.bt.backwardfan[(Int(j),J)])
end
=#

#=
function trace_backwards(s::State, bt::BoundaryTriangulation, c::Candidate)
	i,J = bt.backward[s.e]
	fx = c(inv(J), i+s.x)
	j=Int(floor(fx))
	return State(fx-Int(floor(fx)), bt.forwardfan[(j,J)])
end
=#


function random_cand(bt, thickness, roundmode)
    d=Dict{Junction, DiscreteHomeo{Tuple{Int,Int}}}()
    for j in bt.junctions
        d[j] = random_discrete_homeo(j, thickness, roundmode)
        d[inv(j)] = inv(d[j])
    end

    return Cand(bt, ArrayDict(d))
end


#=
function random_discrete_homeo(j::Junction, thickness::Int, roundmode::RoundMode)
    left_heights = Rational{Int}[i//(thickness * j.left_len) for i in 0:j.left_len * thickness]
    right_heights = Rational{Int}[i//(thickness * j.right_len) for i in 0:j.right_len * thickness]

    ordering = random_ordering(left_heights, right_heights)
    return DiscreteHomeo(left_heights, right_heights, ordering, LEFT, roundmode)
end
=#
function rung_percentage(c::Cand; time=500)
    return [rung_percentage(c,i; time=time) for i in 1:c.bt.ncusps]
end


function slope(c::Cand; time=200)
	return [slope(c, i; time=time) for i in 1:c.bt.ncusps]
end
function uncertain_slope(c::Cand; time=200)
    #return Measurement{Float64}[measurement(exact_slope(c, i),0) for i in 1:c.bt.ncusps]
	return Measurement{Float64}[uncertain_slope(c, i; time=time) for i in 1:c.bt.ncusps]
end
function exact_slope(c::Cand)
    return [exact_slope(c, i) for i in 1:c.bt.ncusps]
end

#workhorse slope function. Just tries to start tracing from a state
function _slope(c::Cand, s::State; time=200)
	#weight=T[0,0]
    W = c.bt.weights
    
	w1=0
	w2=0

	#weight .+= c.bt.weights[s.e]
	x=W[s.e]
	w1+=x[1]
	w2+=x[2]
	while abs(w1) + abs(w2) < time
		s=trace_forwards(s, c)
		#weight .+= c.bt.weights[s.e]
		x=W[s.e]
		w1+=x[1]
		w2+=x[2]
		#println(weight)
	end

	#=
	for i in 1:20
		println(s)
		s=trace_backwards(s, c.bt, c)
		weight = weight + c.bt.weights[s.e]
		println(weight)
	end
	=#

	return (w1,w2)
end

function _exact_slope(c::Cand, s::State{T}) where T
    #return _slope(c, s, time=1000)
    visited = Dict{State{T},Slope}()#stores sum of weights up to, but not including this edge.
    W = c.bt.weights

    w = Slope([0,0])

    while !(haskey(visited,s))
        visited[s] = w
        w += W[s.e]
		s=trace_forwards(s, c)
	end
    return w - visited[s]
end

function rung_percentage(c::Cand{H}, j::Int; time=500) where {H}
    rungs = Set(Iterators.flatten(c.bt.rungs[j]))
    
    sequence = Bool[]
    s=State(rand_init(H), c.bt.rungs[j][1][1])
    for i in 1:time
        s=trace_forwards(s,c)
        push!(sequence, s.e in rungs)
    end
    return sum(sequence) / length(sequence)
end

function show_trace(c::Cand{H}; time=50) where {H}
    for j in c.bt.ncusps
        s=State(rand_init(H), c.bt.rungs[j][1][1])
        @show s
        println("cusp $(j)")
        for i in 1:time
            s=trace_forwards(s, c)
            @show s
        end
    end
end

rand_init(::Type{H}) where {H<:DiscreteHomeo2} = 0//1 
rand_init(c::Cand) = rand_init(typeof(c.d[c.bt.junctions[1]]))
rand_init(::Type{DiscreteHomeo{T}}) where {T} = 1
#rand_init(::Type{CompositeHomeo{H}}) where {H} = rand_init(H)
#rand_init(::Type)=0.31432

function slope(c::Cand{H}, i::Int; time=200) where {H}
    w1,w2 = _slope(c, State(rand_init(c),  c.bt.rungs[i][1][1]), time=time)
	if abs(w1)<=1
		#return NaN
        return sign(w2)*CLIP
	end
	#@show w1,w2
	return clip(w2/w1, CLIP)
end

function uncertain_slope(c::Cand{H}, i::Int; time=200) where {H}
    w1,w2 = _slope(c, State(rand_init(c), c.bt.rungs[i][1][1]), time=time)
	if abs(w1)<=1
		return NaN
	end
	return clip(measurement(w2,1)/w1, CLIP)
end

function exact_slope(c::Cand{H}, i::Int) where {H}
    w1, w2 = _exact_slope(c, State(rand_init(H),  c.bt.rungs[i][1][1]))
	if w1==0
		#return NaN
        return sign(w2)*CLIP
	end
	#@show w1,w2
	return clip(w2/w1, CLIP)
end

function clip(x, bounds)
	ret = max(min(x, bounds), -bounds)
	@assert !isnan(ret)
	@assert !isinf(ret)
	return ret
end


function hasnan(s)
	ret = any(map(isnan,s))
	if ret
		#println("rejecting $(s)")
	end
	return ret
end

struct Longitude
	bt::BoundaryTriangulation
	weights::OffsetArray
end

function push!(e::Envelope, x::Cand)
	push!(e, (exact_slope(x), x))
end


#=
function longitude_to_candidate(bt::BoundaryTriangulation, l)
    d=Dict{Junction, DiscreteHomeo}()
    for J in bt.junctions
        #print_junction(bt, J)
        corresponding_left_heights = Rational{Int}[]
        corresponding_right_heights = Rational{Int}[]

        for k in 0:J.left_len-1
            tri_num, index = bt.forwardfan[(k,J)] :: Track
            for j in 1:l[tri_num]
                push!(corresponding_left_heights, (k + j//(l[tri_num]+1))// J.left_len)
            end
        end


        for k in 0:J.right_len-1
            tri_num, index = bt.backwardfan[(k,J)] :: Track
            for j in 1:l[tri_num]
                push!(corresponding_right_heights, (k + j//(l[tri_num]+1))// J.right_len)
            end
        end

        @assert length(corresponding_left_heights) == length(corresponding_right_heights)

        all_left_heights = sort(vcat(corresponding_left_heights, [i//J.left_len for i in 0:J.left_len]))
        all_right_heights = sort(vcat(corresponding_right_heights, [i//J.right_len for i in 0:J.right_len]))

        push!(corresponding_left_heights, 1//1)
        push!(corresponding_right_heights, 1//1)

        push!(corresponding_left_heights, 0//1)
        push!(corresponding_right_heights, 0//1)

        left_heights = Rational{Int}[]
        right_heights = Rational{Int}[]
        ordering = Dir[]
        while length(all_left_heights) > 0 || length(all_right_heights) > 0
            if !(all_left_heights[1] in corresponding_left_heights)
                push!(left_heights, popfirst!(all_left_heights))
                push!(ordering, LEFT)
            elseif !(all_right_heights[1] in corresponding_right_heights)
                 push!(right_heights, popfirst!(all_right_heights))
                 push!(ordering, RIGHT)
            else
                push!(left_heights, popfirst!(all_left_heights))
                push!(right_heights, popfirst!(all_right_heights))
                push!(ordering, BOTH)
            end
        end
        #=
        @show left_heights, right_heights, ordering
        @show J.left_len
        @show J.right_len
        =#

        H = DiscreteHomeo(left_heights, right_heights, ordering, LEFT, DOWN)
        d[J] = H
        d[inv(J)] = inv(d[J])
    end

    return Cand(bt, ArrayDict(d))
end
=#


function slopes(l::Longitude)
	tmp = [sum(l.weights[i]*l.bt.weights[(i,j)] for (i,j) in edgelist) for edgelist in l.bt.alledges]
	for (ind,i) in enumerate(tmp)
		#=
		for x in [i=> l.weights[i] for i in 0:length(l.weights)-1]
			@show x
		end
		@show i
		@show abs.(i .- round.(Int,i))
		=#
		#c=longitude_to_candidate(bt,l.weights)
		#=
		if sum(abs.(i .- round.(Int,i))) > 0.0001
			@show slope(c, s=State(0.5, bt.firstrungs[1]), verbose=true)
		end
		=#

		#@assert sum(abs.(i .- round.(Int,i))) <= 0.0001
	end
	return [sum(l.weights[i]*l.bt.weights[(i,j)] for (i,j) in edgelist) for edgelist in l.bt.alledges]
end

function relu(x)
	max(x,0)
end


function objective(::Type{Upper}, slope, old_slope, target)
    if any(map(isnan, slope)) || any(map(isinf, slope))
		return -100.0
	end
	@assert length(slope)==length(old_slope)==length(target)
	tmp = sum(relu.(slope.-old_slope) .- 10 * relu.(slope .- target) .- 10 * relu.(old_slope .- slope))
	#tmp = sum((slope.-old_slope) .- 10 * relu.(slope .- target))
	@assert !isnan(tmp)
	@assert !isinf(tmp)
	return tmp
	#return sum(slope)
end
function objective(::Type{Lower}, slope, old_slope, target)
    if any(map(isnan, slope)) || any(map(isinf, slope))
		return -100.0
	end
	@assert length(slope)==length(old_slope)==length(target)
	tmp = sum(relu.(old_slope.-slope) .- 10 * relu.(target .- slope) .- 10 * relu.(slope.-old_slope))
	@assert !isnan(tmp)
	@assert !isinf(tmp)
	return tmp
	#return sum(slope)
end

#=

function Mannealing(f, initial, _jiggle, betastart, betafinish, nsteps; verbose=true, minacc = 100, maxacc = 5000, radius=0.001)
	#linear annealing on range betastart, betafinish
    current = to_mutable(initial)
	curracc = minacc
	currval = f(current; acc=curracc)
	reject_count = 0
	accept_count = 0
    working = to_mutable(initial) #lazy way to copy

	vals=[]
	push!(vals,currval)
	for (i,currbeta) in zip(1:nsteps, range(betastart,betafinish, nsteps))
        copy_to!(current, working)
        jiggle!(working, radius)
		jig = working
		newacc = minacc
		newval = f(jig; acc=minacc)

		r=rand()

		@label here
		dE = exp(currbeta * (-newval+ currval))
		prob = 1/(1+dE)
		#prob will be an interval
	
		#@show prob
		if r > prob
			reject_count += 1
		elseif r < prob
            copy_to!(jig, current)#current = jig
			currval = newval
			curracc = newacc
			accept_count += 1
		elseif (curracc >= maxacc && newacc >= maxacc)
			println("not enough accuracy")
			@show prob, dE
			@show currval
			@show newval
	#		@assert false
		else #improve the accuracy of our computation
			if newacc < maxacc
				newacc *= 3
				newval = f(jig; acc=newacc)
			end
			if curracc < newacc
				curracc *= 3
				currval = f(current; acc = curracc)
			end
			@goto here
		end
		if verbose
			push!(vals, currval)
		end

		if verbose && i%10000 == 0
			@show (exact_slope(current))
			@show (accept_count, reject_count) 
			p=plot(1:length(vals), vals)
			display(p)
		end
	end
	@show (accept_count, reject_count)
    return to_static(current)
end
=#

function annealing(f, initial, jiggle, betastart, betafinish, nsteps; verbose=true, minacc = 100, maxacc = 5000)
	#linear annealing on range betastart, betafinish
	current = initial
	curracc = minacc
	currval = f(current; acc=curracc)
	reject_count = 0
	accept_count = 0

	vals=[]
	push!(vals,currval)
	for (i,currbeta) in zip(1:nsteps, range(betastart,betafinish, nsteps))
		jig = jiggle(current)
		newacc = minacc
		newval = f(jig; acc=minacc)

		r=rand()

		@label here
		dE = exp(currbeta * (-newval+ currval))
		prob = 1/(1+dE)
		#prob will be an interval
	
		#@show prob
		if r > prob
			reject_count += 1
		elseif r < prob
			current = jig
			currval = newval
			curracc = newacc
			accept_count += 1
		elseif (curracc >= maxacc && newacc >= maxacc)
			println("not enough accuracy")
			@show prob, dE
			@show currval
			@show newval
	#		@assert false
		else #improve the accuracy of our computation
			if newacc < maxacc
				newacc *= 3
				newval = f(jig; acc=newacc)
			end
			if curracc < newacc
				curracc *= 3
				currval = f(current; acc = curracc)
			end
			@goto here
		end
		if verbose
			push!(vals, currval)
		end

		if verbose && i%10000 == 0
			@show (exact_slope(current))
			@show (accept_count, reject_count) 
			p=plot(1:length(vals), vals)
			display(p)
		end
	end
	@show (accept_count, reject_count)
	return current
end

function try_improve(E::Envelope{S,T,D}; nsubdivide=0, iters=50000, time=1000, target = [1000,1000], beta=500, radius=0.001, min_cands=3) where {S,T,D}
	accurate_E=Envelope{S,T,D}()
    
    cands = sort([(objective(S, x[1], x[1], target), x[2]) for x in E.A], by=x->-x[1])

    maxind = min(min_cands, length(cands))
    while maxind < length(cands) && cands[maxind+1][1] > -0.001
        maxind += 1
    end

    println("restricted to $maxind / $(length(cands))")

	@threads for i in 1:maxind
		_, oldcand = cands[i]
		old_v = exact_slope(oldcand)

		for j in 1:nsubdivide
			oldcand = subdivide(oldcand)
		end

		push!(accurate_E, (old_v, oldcand))
        #=
		if objective(S, old_v, old_v, target) < -0.001
			continue
		end
        =#

		@time newcand = annealing((c; acc=100)->objective(S, uncertain_slope(c,time=acc), old_v, target), oldcand, c->jiggle(c,radius), beta, beta, iters; verbose=false)
		new_v = exact_slope(newcand)
        #@show old_v, new_v 
        
        lock(stdout) 
        begin
            print(old_v)
            print("->[")
            for (x,val) in zip(sign.(new_v .- old_v), new_v)
                a = S==Upper ? 1 : -1
                if x * a > 0
                    col = :green
                elseif x * a < 0
                    col = :red
                else
                    col = :yellow
                end
                printstyled(val, color = col)
                print(",")
            end
            print("]")
            println()
        end
        unlock(stdout)
		push!(accurate_E, (new_v, newcand))
	end
	return accurate_E
end


function random_trials(bt; thickness=8, roundmode=DOWN, ntrials=100000)
	E=PEnvelope()
	#trials=[T[] for i in 1:N]


	#inrange(x) = all(-5 < i < 5 for i in x)

	@threads for i in 1:ntrials
		c=random_cand(bt,thickness,roundmode)
		#trials[i] = slope(c)
		push!(E, (exact_slope(c),c))
	end

	#@show [x[1] for x in E.A]
	#vals = collect(filter(inrange, trials))

	#accurate_E = try_improve(E)
	#accurate_envelope = sort(collect(filter(inrange, [v for (v,c) in accurate_E.A])))

	#p=scatter([x[1] for x in vals], [x[2] for x in vals], markersize=2)
	#scatter!(p, [x[1] for x in accurate_envelope], [x[2] for x in accurate_envelope])
	#display(p)
	return E
end

function extreme_candidates(bt)
	Elower=Envelope{Lower}()
	Eupper=Envelope{Upper}()

	@threads :greedy for tmp in Iterators.product(([Piecewise(Linear(),Linear(),0.01, 0.99), Piecewise(Linear(),Linear(),0.99, 0.01), Linear()] for i in 1:length(bt.junctions))...)
		c=Candidate(bt, Dict{Junction, Union{Linear,Piecewise}}(j=>f for (j,f) in zip(bt.junctions,tmp)))
		ss = exact_slope(c)
		push!(Elower, (ss,c))
		push!(Eupper, (ss,c))
	end
	return Elower, Eupper
end

#=
function scatter_envelope!(p, E::Envelope)
	vals = [x[1] for x in E.A]
	@assert length(vals) > 0
	scatter!(p, [[x[i] for x in vals] for i in 1:length(vals[1])]..., markersize=2)
end
=#


function inbounds(pt)
	return all(abs.(pt) .<= CLIP)
end

function staircase(E::Envelope{Upper})
	pts = [x[1] for x in E.A]
	pts=filter(inbounds, sort!(pts, by=x->x[1]))
	return sort(vcat(pts,[(pts[i][1], pts[i+1][2]) for i in 1:length(pts)-1]), by=x->(x[1],-x[2]))
end

function crevices(E::Envelope{Upper})
	pts = [x[1] for x in E.A]
	pts=filter(inbounds, sort!(pts, by=x->x[1]))
    return [[pts[i][1], pts[i+1][2]] for i in 1:length(pts)-1]
end

function staircase(E::Envelope{Lower})
	pts = [x[1] for x in E.A]
	pts = filter(inbounds, sort!(pts, by=x->x[1]))
	return sort(vcat(pts,[(pts[i+1][1], pts[i][2]) for i in 1:length(pts)-1]), by=x->(x[1],-x[2]))
end

function crevices(E::Envelope{Lower})
	pts = [x[1] for x in E.A]
	pts = filter(inbounds, sort!(pts, by=x->x[1]))
    return [[pts[i+1][1], pts[i][2]] for i in 1:length(pts)-1]
end


function normalizedchi(L::Longitude)
	ss=slopes(L)
	chi = -sum(L.weights)//2
	npunctures = [gcd(a,b) for (a,b) in ss]
	g=gcd(npunctures...)
	closed_chi = chi + sum(npunctures)#npunctures[2] + npunctures[1]
	#return closed_chi/g #this is a bit suspicious. How do we know that this is the multiplicity of the surface?
    return closed_chi
end

function updateith(A::Vector{T}, i::Int, val::T) where {T}
    B = copy(A)
    B[i]=val
    return B
end

function constraints_multi(L::Longitude)
    ss=slopes(L)
    chi = Int(-sum(L.weights)//2) #Euler characteristic of the punctured surface

    npunctures = [gcd(a,b) for (a,b) in ss]
end

function constraints(L::Longitude)
	ss = slopes(L)
    rats = map(slope_to_rat, ss)
    ncusps = length(ss)
    chi = Int(-sum(L.weights)//2) #Euler characteristic of the punctured surface
	#@show chi


	npunctures = [gcd(a,b) for (a,b) in ss]
	#@show npunctures
    closed_chi = chi + sum(npunctures)

    @show chi
    @show npunctures

	#=
	if closed_chi >= 0
		#@show L
		@show ss
		@show closed_chi
	end
	=#


    #=
	q,p = ss[1]
	s,r = ss[2]

    @assert q>=0
    @assert s>=0
    =#

    nladders = map(x->length(x)//2, L.bt.rungs)

    #prongcounts = (q*nladders[1]//gcd(p,q), s*nladders[2]//gcd(r,s))

    prongcounts = [ss[i][1] * nladders[i] //gcd(ss[i][1], ss[i][2]) for i in 1:ncusps]

    @assert closed_chi == -sum(npunctures[i]*(prongcounts[i]-2) for i in 1:ncusps)//2

    #@assert closed_chi==-sum(npunctures[1]*(prongcounts[1]-2) + npunctures[2]*(prongcounts[2]-2))//2
	#=
		if q-closed_chi == 0 
			return (0,0)
		else
			return (p/(q - closed_chi),  r/s)
		end
	=#

	#@show (q,p,s,r)
	#@show npunctures
	ret = []
	#suppose (0,0) is an S^1 \times S^2 surgery 
	#
	
	#@show chi+npunctures[1]
	#@show chi+npunctures[2]
    #
    #

    for i in 1:ncusps
        bounds, info = constraint3(ss[i], chi + sum(npunctures) - npunctures[i])
        for (j,bound) in enumerate(bounds)
            info_full = (info...,npunc=npunctures[i], int_prongs = [prongcounts[i] for j in 1:ncusps if i !=j], ext_prongs=prongcounts[i], longitude=L.weights, ss=ss, dir=(i, j==1 ? -1 : 1))
            push!(ret, (updateith(rats, i, slope_to_rat(bound)), info_full))
        end
    end

    #=

    #if npunctures[2]==1
        x,info = constraint3(r,s,chi+npunctures[1], upper=true)
        push!(ret, ((p//q, x), (info...,npunctures=npunctures[2],interior_prong = prongcounts[1], exteriorprong=prongcounts[2], longitude=L.weights, ss=ss, dir=:up)))

        x,info = constraint3(r,s,chi+npunctures[1], upper=false)
        push!(ret, ((p//q, x),(info...,npunctures=npunctures[2],interior_prong = q*nladders[1]//gcd(p,q), exteriorprong=prongcounts[2],longitude=L.weights, ss=ss, dir=:down)))
    #end
    #if npunctures[1]==1
        x,info = constraint3(p,q,chi+npunctures[2], upper=true)
        push!(ret, ((x, r//s),(info...,npunctures=npunctures[1],interior_prong=s*nladders[2]//gcd(r,s), exteriorprong=prongcounts[1],lonigtude=L.weights, ss=ss, dir=:right)))

        x,info = constraint3(p,q,chi+npunctures[2], upper=false)
        push!(ret, ((x, r//s), (info...,npunctures=npunctures[1],interior_prong=s*nladders[2]//gcd(r,s), exteriorprong=prongcounts[1],longitude=L.weights,ss=ss, dir=:left)))
    #end
    =#

    #=
    if 0 in Iterators.flatten(ret)
        println("----")
        @show npunctures
        @show ss
        @show (p//q,constraint3(r,s,chi+npunctures[1]; verbose=true))
        @show (p//q, -constraint3(-r,s,chi+npunctures[1]; verbose=true))
        @show (constraint3(p,q,chi+npunctures[2]; verbose=true), r//s)
        @show (-constraint3(-p,q,chi+npunctures[2]; verbose=true), r//s)
    end
    =#

	return ret

	#=
	for sgn in [-1,1]
		if s - (chi+npunctures[1]) != 0

			n=floor(r/s)
			rp,sp = numerator(r/s-n), denominator(r/s-n)

			push!(ret,(p/q,  n + r/(s - (chi+npunctures[1]))))

		else
			push!(ret,(CLIP,CLIP))
		end
		if q + sgn*(chi+npunctures[2]) != 0
			push!(ret, (p/(q + sgn*(chi+npunctures[2])), r/s))
		else
			push!(ret,(CLIP,CLIP))
		end
	end
	return ret

	=#

	#r/s, r/(s+chi)

	#(r+s)/s, (r+s)/(s+chi)
end

function intersection_weights(bt::BoundaryTriangulation, loop::Vector{Track})
	weights=DefaultDict{Track,Int}(0)
	for i in 1:length(loop)
		#@show loop[i]
		t1,t2 = loop[i], loop[mod1(i+1, length(loop))]
		i1, J1=bt.forward[t1]
		i2, J2=bt.backward[t2]
		@assert J1==J2
		J=J1

		#@show i1, i2

		for k in 0:J.left_len-1
			#@show (k,J)
			#@show bt.forwardfan[(k,J)]
			weights[bt.forwardfan[(k, J)]]+= if k < i1
				1#1//2
			elseif k > i1
				0#-1//2
			else
				0
			end
		end
		for k in 0:J.right_len-1
			#@show (k,J)
			#@show bt.backwardfan[(k,J)]
			weights[bt.backwardfan[(k, J)]] += if k < i2
				-1#-1//2
			elseif k > i2
				0#1//2
			else
				0
			end
		end
	end
	return weights
end

function ladderpoles(bt::BoundaryTriangulation)
    return [ladderpoles(bt,i) for i in 1:bt.ncusps]
end

function ladderpoles(bt::BoundaryTriangulation, i::Int)
    return [ 
    begin
        _,J=bt.backward[rungs[1]]
        f(bt,J)
    end

    for (rungs,f) in zip(bt.rungs[i], Iterators.cycle((upward_ladderpole, downward_ladderpole)))]
end

function upward_ladderpole(bt::BoundaryTriangulation, J::Junction)
    tracks = Track[]
    currJ = J

    while true
        curr_edge = bt.backwardfan[(currJ.right_len-1,currJ)]
        push!(tracks, curr_edge)
        _,currJ = bt.forward[curr_edge]

        if currJ == J
            return tracks
        end
    end
end

function downward_ladderpole(bt::BoundaryTriangulation, J::Junction)
    tracks = Track[]
    currJ = J

    while true
        curr_edge = bt.forwardfan[(currJ.left_len-1,currJ)]
        push!(tracks, curr_edge)
        _,currJ = bt.backward[curr_edge]

        if currJ == J
            return tracks
        end
    end
end


#=
function trace_forwards(s::State{T}, c::Cand) where {T}
    #@show s.x
    #@show s.e
	@assert 0<= s.x <=1
	i,J = c.bt.forward[s.e]
	
#	print_junction(bt,J)
    fx = c(J,i+s.x)
    #@show fx


    j=Int(floor(fx))
	#@show j
	
	if j==J.right_len
        @assert fx==j
		#fx = j
		j=j-1
	end

    #println("next\n\n")

    return State{T}(fx-j, c.bt.backwardfan[(Int(j),J)])
end
=#


function compute_loop(bt,rungs) 
    #rungs is a list of rungs, all in the first ladder. Returns a loop in that cusp, intersection number one with that ladder
	S=Set{Track}(rungs)
	Q=Queue{Vector{Track}}()
	for rung in circshift(rungs,-1) #the second rung is the bottommost
		enqueue!(Q, Track[rung])
	end
	while !isempty(Q)
		l=dequeue!(Q)

		_,J=bt.forward[l[end]]
		for k in 0:J.right_len-1
			next = bt.backwardfan[(k, J)]
			if next == l[1]
				return l
			end
			if !(next in S) && !(next in l)
				enqueue!(Q,vcat(l, Track[next]))
			end
		end
	end
	@assert false #failed to find a loop
end

function compute_weights!(bt::BoundaryTriangulation)
    for runglist in bt.rungs
        push!(bt.zero_loops, compute_loop(bt, runglist[1]))
    end

    for l in bt.rungs
        for track in l[1] #only assign weight 1 to the first ladder.
            bt.weights[track] = Slope([1,0])
        end
    end

	for loop in bt.zero_loops
		#@show loop
		weights = intersection_weights(bt, loop)
		for (k,v) in weights
			bt.weights[k] = Slope([bt.weights[k][1], v])
		end
	end


end

function degen_to_snappy_basis_change(bt::BoundaryTriangulation)
    ret = Array{Int,2}[]
    for i in 1:bt.ncusps
        degen = bt.poles[i][1]
        zero = bt.zero_loops[i]

        
        m1, l1 = sum(bt.snappy_weights[e] for e in zero)
	    m2, l2 = sum(bt.snappy_weights[e] for e in degen)

        push!(ret, Int[-l1 -l2; m1 m2])
    end
    @assert all(abs(det(A))==1 for A in ret)
    
    return ret
end

function snappy_to_degen_basis_change(bt::BoundaryTriangulation)
    return Array{Int,2}[inv(A) for A in degen_to_snappy_basis_change(bt)]
end

#=
function constraint2(r::Int, s::Int, chi)
	#r/s is the fiber slope
	#the degeneracy slope is 1/0
	
	@assert s>=0
	if s==0
		return NaN
	end
	@assert s>0

	g=gcd(r,s)
	r=r//g
	s=s//g
	chi=chi//g

	if s+chi == 0
		return NaN
	elseif chi>=0
		return r//s
	elseif s+chi > 0
		n=floor(r//s)

		rprime = r-n*s
		sprime = s

		iota = sprime #intersection number between fiber and degeneracy

		m=ceil(chi//iota)
		@show chi,iota,m

		#=
		(0,1)/ + (s,r)/chi

		(0,1)/iota + (s,r)/(m*iota)
		degeneracy/iota + fiber/(m*iota)

		(0,1)*m + (sprime,rprime)
		=#

		return n//1+rprime//(s+m)
		#when chi is big, it should return something close to the degeneracy slope
		#when chi is small, it should return close to the fiber slope
	else 
		return NaN
	end

end
=#


#=
function constraint(r::Int,s::Int,chi) #careful, it's possible gcd(r,s)=/= 1
	#appears that this is correct, at least when r=1 mod s
	
	@assert s>=0
	if s==0
		return NaN
	end
	@assert s > 0
	n=floor(r//s)

	#=
	if abs(chi) > abs(s)/2
		return NaN
	end
	=#

	if s+chi == 0
		return NaN
	elseif chi>=0
		return r//s
	elseif s+chi > 0
		return n//1 + (r-n*s)//(s + chi)
		#return r//(s + chi)
	else
		return NaN
	end
end
=#

function n_orthogonal(v::AbstractVector{Int}) #Find a determinant 1 matrix that sends v to [k,0]
    @assert length(v)==2
    g = gcd(v[1], v[2])
    A = n_orthogonal(v, [1 0; 0 1])
    @assert A*v == [g,0]
    @assert A[1,1]*A[2,2]-A[1,2]*A[2,1]==1
    return A
end

function n_orthogonal(v::AbstractVector{Int},A::Matrix{Int})
    r,s = v
    if v[2]==0 && v[1] > 0
        return A
    else
        if r < 0 && s < 0 
            B=[-1 0; 0 -1]
        elseif r <= 0
            B=[0 1; -1 0]
        elseif s < 0
            B=[0 -1; 1 0]
        else
            @assert r>=0
            @assert s>=0
            if r >= s
                B = [1 -1; 0 1]
            else 
                B = [1 0; -1 1]
            end
        end
        return n_orthogonal(B*v, B*A)
    end
end

function constraint3(fiber_slope::Slope, chi::Int; verbose=false)
    @show fiber_slope, chi

    g=gcd(fiber_slope...,chi)
    fiber_slope = Slope(fiber_slope .// g)
    chi = Int(chi//g)

    #=
    if s==0
        return (NaN, "fractional DT = 0")
    end
    =#
    #=
    if chi >= 0
        return (r//s, "chi>=0")
    end
    =#

    A = n_orthogonal(fiber_slope)
    Ainv = round.(Int,inv(A))

    degen_slope = A * [0,1]
    degen = degen_slope[1]//degen_slope[2]

    if chi > 0
        return [Ainv*[-1,0], Ainv*[1,0]], (chi=chi, degen_slope=degen_slope, upper=1//0, lower=-1//0, A=A)
    end

    upperbound = minimum(Int(ceil(degen*chi0))//chi0 for chi0 in 1:max(-chi,1))
    lowerbound = maximum(Int(floor(degen*chi0))//chi0 for chi0 in 1:max(-chi,1))

    if verbose
        @show A
        @show fiber_slope
        @show chi
        @show degen_slope
        @show degen
        @show upperbound
        @show lowerbound
    end

    return [Ainv*[lowerbound,1],Ainv*[upperbound,1]],  (chi=chi, degen_slope=degen_slope, upper=upperbound, lower=lowerbound, A=A)
end

function bound(E::Envelope{Upper}, r)
	pts = filter(pt-> pt[1] > r+0.0000001, [x[1] for x in E.A])
	return maximum(x[2] for x in pts; init=-CLIP)
end
function bound(E::Envelope{Lower}, r)
	pts = filter(pt-> pt[1] < r-0.000001, [x[1] for x in E.A])
	return minimum(x[2] for x in pts; init=CLIP)
end

function bound(E::Envelope{Upper}, s::Tuple, i::Int) #s is a tuple of slopes, and i is the coordinate in which we want to maximize
    @assert length(s)==2

    if i==2
        pts = filter(pt-> pt[1] > s[1]+0.0000001, [x[1] for x in E.A])
        return maximum(x[2] for x in pts; init=-CLIP)
    else
        @assert i==1
        pts = filter(pt-> pt[2] > s[2]+0.0000001, [x[1] for x in E.A])
        return maximum(x[1] for x in pts; init=-CLIP)
    end
end
function bound(E::Envelope{Lower}, s::Tuple, i::Int) #s is a tuple of slopes, and i is the coordinate in which we want to maximize
    @assert length(s)==2

    if i==2
        pts = filter(pt-> pt[1] < s[1]-0.000001, [x[1] for x in E.A])
        return minimum(x[2] for x in pts; init=CLIP)
    else
        @assert i==1
        pts = filter(pt-> pt[2] < s[2]-0.000001, [x[1] for x in E.A])
        return minimum(x[1] for x in pts; init=CLIP)
    end
end

#=
function components(L::Longitude)
	for J in L.bt.junctions
		forward = []
		backward = []
		for k in 0:J.left_len-1
			for t in L.bt.backwardfan[(J,k)]
			end
		end
	end
end
=#



#Improved setup
#We have an edge height type (indicating the height inside an edge, or more properly, inside a triangle)
#And a junction height type (indicating the height inside the left or right fan)
#And injection maps (going from (edge, height) junction height
