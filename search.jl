const T=Float64
using OffsetArrays
import Base: inv, getindex, setindex!, hash, push!, length
using DataStructures
using Plots
using Profile
#using ProfileView
using Base.Threads
using PlotlyJS
using Measurements
using StaticArrays

#plotlyjs()

import Plots: plot


abstract type Homeo end

#really, the identity function
struct Linear <: Homeo
end

const Track = Tuple{Int,Int} #edge index and vertex number
#represents a piecewise linear function from [0,1] to [0,1]
struct Piecewise <: Homeo
	left::Union{Piecewise, Linear}
	right::Union{Piecewise, Linear}
	x::T
	fx::T #x maps to fx
end

abstract type Comp
end

struct Upper <: Comp
end

struct Lower <: Comp
end

struct Eq <: Comp
end


struct Junction
	index::Int
	inv::Bool
	left_len::Int
	right_len::Int
end


cartesianlength(::Type{Junction}) = 2
cartesianlength(::Type{Tuple{S,T}}) where {S,T} = cartesianlength(S) + cartesianlength(T)
cartesianlength(::Type{Int}) = 1
tupleIndex(J::Junction) = (J.index+1, Int(J.inv)+1)
#tupleIndex((x,y)::Tuple{S,T}) where {S,T} = (tupleIndex(x)..., tupleIndex(y)...)
tupleIndex((x,y)::Tuple{Int,Int}) = (x+1,y+1)
tupleIndex(x::Int) = x+1
tupleIndex((x,J)::Tuple{Int,Junction}) = (x+1, J.index, Int(J.inv)+1)

struct ArrayDict{S,T,N}
	data::Array{T,N}
end

function getindex(A::ArrayDict{S,T,N}, i::S) where {S,T,N}
	@assert isassigned(A.data, tupleIndex(i)...)
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
	A=Array{T,N}(undef, length.(ranges)...)	
	data = OffsetArray(A, ranges...)
	ret = ArrayDict{S,T,N}(A)

	for (k,v) in D
		ret[k]=v
	end
	return  ret
end

struct BoundaryTriangulation
	forward::ArrayDict{Track, Tuple{Int,Junction}, 2}
	backward::ArrayDict{Track, Tuple{Int,Junction}, 2}
	forwardfan::ArrayDict{Tuple{Int,Junction}, Track, 3}
	backwardfan::ArrayDict{Tuple{Int,Junction}, Track, 3}
	weights::ArrayDict{Track, SVector{2,T}, 2}
	junctions::Vector{Junction}
	firstrungs::Vector{Track}#one rung in each cusp
	alledges::Vector{Vector{Track}}#edges appearing in each cusp
	#fans is a list of pairs of fans
	
	function BoundaryTriangulation(fans, face_coorientations, firstrungs, alledges, _weights)
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

		weights=Dict{Track,SVector{2,T}}()
		for (x,y) in _weights
			weights[x] = y
		end
		return new(ArrayDict(forward),ArrayDict(backward),ArrayDict(forwardfan),ArrayDict(backwardfan),ArrayDict(weights),junctions,firstrungs, alledges)
	end
end

struct Candidate
	bt::BoundaryTriangulation
	d::Dict{Junction, Union{Piecewise,Linear}}
end
struct Envelope{S} #keep track of local maxes
	A::Vector{Tuple{Vector{T},Candidate}}
	L::SpinLock
end

function prunings(p::Piecewise)
	chnl = Channel{Homeo}(3)
	put!(chnl, Linear())
	@async begin
		for i1 in prunings(p.left)
			push!(chnl, Piecewise(i1,p.right,p.x,p.fx))
		end
		for i2 in prunings(p.right)
			push!(chnl, Piecewise(p.left,i2,p.x,p.fx))
		end
		close(chnl)
	end
	return chnl
end

function complexity(p::Piecewise)
	return complexity(p.left) + complexity(p.right)
end
function complexity(p::Linear)
	return 1
end

function prunings(l::Linear)
	return (Linear(),)
end


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

function inv(p::Piecewise)
	Piecewise(inv(p.left),inv(p.right),p.fx,p.x)
end

function inv(x::Linear)
	Linear()
end

function (f::Piecewise)(y::T)
	@assert 0 <= y <=1
	if y < f.x
		ret = f.left(y/f.x) * f.fx
	else
		ret = f.right((y-f.x)/(1-f.x)) * (1-f.fx) + f.fx
	end
	@assert 0 <= ret <= 1
	return ret
end

#=
function (f::Optimistic)(y::T)
	for (h,b) in f.order
		if !b && 

end
=#


function (f::Linear)(y::T)
	@assert 0 <= y <= 1
	return y
end

function cwise((fnum,edge_index)::Tuple{Int,Int}, coor::Int)
	(fnum,mod(edge_index-coor,3))
end
function ccwise((fnum,edge_index)::Tuple{Int,Int}, coor::Int)
	(fnum,mod(edge_index+coor,3))
end


function Base.hash(J::Junction, h::UInt)
	return hash(J.index, hash(J.inv, h))
end


struct State
	x::T#height in the edge
	e::Track#which track we're on
end

function inv(e::Junction)
	return Junction(e.index,!e.inv,e.right_len,e.left_len)
end



function prunings(c::Candidate)
	ch = Channel{Candidate}(3)
	D=Dict{Junction, Union{Piecewise,Linear}}()
	for J in c.bt.junctions
		D[J] = c[J]
	end
	@async begin
		for J in c.bt.junctions
			for p in prunings(D[J])
				if complexity(p) < complexity(D[J])
					cnew=Candidate(c.bt,Dict{Junction,Union{Piecewise,Linear}}())
					for _J in bt.junctions
						cnew[_J]=c[_J]
					end
					cnew[J]=p
					put!(ch, cnew)
				end
			end
		end
		close(ch)
	end
	return ch
end

function complexity(c::Candidate)
	return sum(complexity(c[J]) for J in c.bt.junctions)
end

function prune(c::Candidate)
	val = approximant_all_slopes(c)
	@show complexity(c)

	@label here
	for cnew in prunings(c)
		if approximant_all_slopes(cnew,time=10000) == val
			#@show complexity(cnew)
			c=cnew
			@goto here
		end
	end

	@show complexity(c)
	println("--")

	return c
end

function prune(E::Envelope{S}) where {S}
	return Envelope{S}([(val, prune(c)) for (val,c) in E.A])
end

function prune!(E::Envelope)
	@threads for i in eachindex(E.A)
		val,c = E.A[i]
		E.A[i] = (val, prune(c))
	end
end

function getindex(c::Candidate, J::Junction)
	if J.inv && !haskey(c.d, J)
		c.d[J] = inv(c[inv(J)])
	end
	return c.d[J]
end

function setindex!(c::Candidate, f::Union{Piecewise,Linear}, J::Junction)
	c.d[J] = f
	c.d[inv(J)] = inv(f)
end

function (c::Candidate)(J::Junction,x)
	@assert 0 <= x <= J.left_len
	return c[J](x/J.left_len) * J.right_len
end

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

function trace_forwards(s::State, bt::BoundaryTriangulation, c::Candidate)
	@assert 0<= s.x <=1
	i,J = bt.forward[s.e]
	
	#print_junction(bt,J)
	fx = c(J,i+s.x)


	j=Int(floor(fx))
	#@show j
	
	#=
	if !haskey(junction_crossings, J)
		junction_crossings[J]=Set()
	end
	push!(junction_crossings[J], (i,j))
	=#
	if j==J.right_len
		@assert fx - j < 0.00000001
		j=j-1
	end
	#=
	if !haskey(bt.backwardfan, (j,J))
		@show s
		@show j
		@show fx
		@show keys(bt.backwardfan)
		@show J.right_len
	end
	=#

	return State(fx-Int(floor(fx)), bt.backwardfan[(j,J)])
end

function trace_backwards(s::State, bt::BoundaryTriangulation, c::Candidate)
	i,J = bt.backward[s.e]
	fx = c(inv(J), i+s.x)
	j=Int(floor(fx))
	return State(fx-Int(floor(fx)), bt.forwardfan[(j,J)])
end

function random_piecewise(depth)
	if depth<=0
		Linear()
	else
		Piecewise(random_piecewise(depth-1), random_piecewise(depth-1), uniform(0.01,0.99), uniform(0.01,0.99))
	end
end

function subdivide(p::Piecewise)
	return Piecewise(subdivide(p.left), subdivide(p.right), p.x, p.fx)
end
function subdivide(p::Linear)
	return Piecewise(Linear(), Linear(), 0.5, 0.5)
end
function subdivide(c::Candidate)
	return Candidate(c.bt, Dict{Junction, Union{Piecewise,Linear}}(x=>subdivide(y) for (x,y) in c.d))
end

function piecewise(inputs::Vector, outputs::Vector)
	@assert length(inputs)==length(outputs)
	for x in inputs
		@assert 0 <= x <= 1
	end
	for x in outputs 
		@assert 0 <= x <= 1
	end
	@assert inputs == sort(inputs)
	@assert outputs == sort(outputs)

	if length(inputs)==0
		return Linear()
	else
		i = round(Int,ceil(length(inputs)/2))
		ret = Piecewise(piecewise(inputs[1:i-1] ./ inputs[i], outputs[1:i-1] ./ outputs[i]),
						 piecewise((inputs[i+1:end] .- inputs[i]) ./ (1-inputs[i]), (outputs[i+1:end] .- outputs[i]) ./ (1-outputs[i])),
						 inputs[i],
						 outputs[i])
		for (i,j) in zip(inputs, outputs)
			@assert abs(j-ret(i))<0.0001
		end
		return ret
	end
end

function random_candidate(bt,depth)
	c=Candidate(bt, Dict{Junction, Union{Linear,Piecewise}}())
	for j in bt.junctions
		#@show depth + round(Int, log2(j.left_len + j.right_len))
		c.d[j]=random_piecewise(depth + round(Int, log2(j.left_len + j.right_len)))
	end
	return c
end

function uniform(x,y)
	@assert y>=x
	return rand()*(y-x) + x
end

function approximant_all_slopes(c::Candidate; time=10000)
	return [approximant(x,time) for x in all_slopes(c; time=time)]
end

function approximant(x, time)
	if x==Inf
		return x
	elseif x==Inf || x==-Inf
		return 1//0
	elseif x < 0
		-approximant(-x, time)
	elseif abs(x) > 1
		1//approximant(1/x, time)
	else
		rationalize(x; tol=10/time)
	end
end

function all_slopes(c::Candidate; time=200)
	return [slope(c; s=State(0.023423,rung), time=time) for rung in c.bt.firstrungs]
end

function all_uncertain_slopes(c::Candidate; time=200)
	return [uncertain_slope(c; s=State(0.023423,rung), time=time) for rung in c.bt.firstrungs]
end

#=
function all_LazySlopes(c::Candidate)
	return [LazySlope(c, rung) for rung in c.bt.firstrungs]
end


mutable struct LazySlope
	c::Candidate
	weight::Vector{T}
	s::State
end

function value(l::LazySlope)
	return interval(l.weight[2]-1, l.weight[2]+1)/l.weight[1]
end

function LazySlope(c::Candidate, s=State(0.312423, (1,0)))
	return LazySlope(c, T[0,0], s)
end

function compute_slope!(l::LazySlope; time=100)
	while abs(l.weight[1]) + abs(l.weight[2]) < time
		l.s=trace_forwards(l.s, c.bt, c)
		l.weight .+= c.bt.weights[l.s.e]
	end
end
=#


#todo: return confidence interval, and allow to improve the confidence interval with more work.
function slope(c::Candidate; time=200, s=State(0.312423, (1,0)))
	#weight=T[0,0]
	w1=T(0)
	w2=T(0)

	#weight .+= c.bt.weights[s.e]
	x=c.bt.weights[s.e]
	w1+=x[1]
	w2+=x[2]
	while abs(w1) + abs(w2) < time
		#println(s)	
		s=trace_forwards(s, c.bt, c)
		#weight .+= c.bt.weights[s.e]
		x=c.bt.weights[s.e]
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
	return clip(w2/w1, 10)
end

function uncertain_slope(c::Candidate; time=200, s=State(0.312423, (1,0)))
	#weight=T[0,0]
	w1=T(0)
	w2=T(0)
	#weight .+= c.bt.weights[s.e]
	#somehow the above allocates, so we'll be more explicit
	x=c.bt.weights[s.e]
	w1+=x[1]
	w2+=x[2]
	while abs(w1) + abs(w2) < time
		#println(s)	
		s=trace_forwards(s, c.bt, c)
		#weight .+= c.bt.weights[s.e]
		x,y=c.bt.weights[s.e]
		w1+=x
		w2+=y
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
	return clip(measurement(w2,1)/w1, 10)
end

function clip(x, bounds)
	ret = max(min(x, bounds), -bounds)
	@assert !isnan(ret)
	@assert !isinf(ret)
	return ret
end

function Envelope{S}() where {S <: Comp}
	return Envelope{S}(Tuple{Vector{T},Candidate}[], SpinLock())
end

function Envelope{S}(A::Vector) where {S<: Comp}
	return Envelope{S}(A, SpinLock())
end

function Envelope()
	return Envelope{Upper}()
end

function PEnvelope()
	return Envelope{Eq}()
end

function strict_compare(x::Vector{T},y::Vector{T})
	return all(x .<= y)
end

function comp(S::Type{Upper}, x, y)
	return strict_compare(x, y)
end
function comp(S::Type{Lower}, x, y)
	return strict_compare(y, x)
end
function comp(S::Type{Eq}, x, y)
	return x==y
end

function push!(e::Envelope{S}, x::Tuple{Vector{T},Candidate}) where {S}
	lock(e.L) do
		if !any(comp(S, x[1], y[1]) for y in e.A)
			filter!(y->!comp(S,y[1],x[1]), e.A)
			push!(e.A, x)
		end
	end
end

function push!(e::Envelope{S}, x::Tuple{Union{NTuple{N,R}, Vector{R}},Candidate}) where {S, N, R <: Real}
	push!(e, (T[x[1]...], x[2]))
end

function push!(e::Envelope, x::Candidate)
	push!(e, (approximant_all_slopes(x), x))
end

struct Longitude
	bt::BoundaryTriangulation
	weights::OffsetArray
end

function slopes(l::Longitude)
	[round.(Int,sum(l.weights[i]*bt.weights[(i,j)] for (i,j) in edgelist)) for edgelist in bt.alledges]
end

function relu(x)
	max(x,0)
end


function objective(::Type{Upper}, slope, old_slope, target)
	tmp = sum(relu.(slope.-old_slope) .- 10 * relu.(slope .- target) .- 10 * relu.(old_slope .- slope))
	@assert !isnan(tmp)
	@assert !isinf(tmp)
	return tmp
	#return sum(slope)
end
function objective(::Type{Lower}, slope, old_slope, target)
	tmp = sum(relu.(old_slope.-slope) .- 10 * relu.(target .- slope) .- 10 * relu.(slope.-old_slope))
	@assert !isnan(tmp)
	@assert !isinf(tmp)
	return tmp
	#return sum(slope)
end

function parallel_try_improve(E::Envelope; nsubdivide=0, iters=50000, time=1000, target=[1000,1000], beta=500)

end

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
				newacc *= 4
				newval = f(jig; acc=newacc)
			end
			if curracc < newacc
				curracc *= 4
				currval = f(current; acc = curracc)
			end
			@goto here
		end
		push!(vals, currval)

		if verbose && i%10000 == 0
			@show (approximant_all_slopes(current))
			@show (accept_count, reject_count) 
			p=plot(1:length(vals), vals)
			display(p)
		end

	end
	@show (accept_count, reject_count)
	return current
end

function try_improve(E::Envelope{S}; nsubdivide=0, iters=50000, time=1000, target = [1000,1000], beta=500, radius=0.001) where {S}
	@show length(E.A)
	accurate_E=Envelope{S}()
	@threads for i in 1:length(E.A)
		_,oldcand = E.A[i]
		old_v = T.(approximant_all_slopes(oldcand))

		for i in 1:nsubdivide
			oldcand = subdivide(oldcand)
		end

		push!(accurate_E, (old_v, oldcand))
		if objective(S, old_v, old_v, target) < -0.001
			continue
		end

		newcand = annealing((c; acc=100)->objective(S, all_uncertain_slopes(c,time=acc), old_v, target), oldcand, c->jiggle(c,radius), beta, beta, iters; verbose=false)
		new_v = T.(approximant_all_slopes(newcand))
		@show old_v, new_v
		push!(accurate_E, (new_v, newcand))
	end
	return accurate_E
end

length(e::Envelope) = length(e.A)

function random_trials(bt; nsubdivide=2)
	N=100000
	E=PEnvelope()
	trials=[T[] for i in 1:N]


	#inrange(x) = all(-5 < i < 5 for i in x)

	@threads for i in 1:N
		c=random_candidate(bt,nsubdivide)
		trials[i] = all_slopes(c)
		push!(E, (trials[i],c))
	end

	@show [x[1] for x in E.A]
	#vals = collect(filter(inrange, trials))

	#accurate_E = try_improve(E)
	#accurate_envelope = sort(collect(filter(inrange, [v for (v,c) in accurate_E.A])))

	#p=scatter([x[1] for x in vals], [x[2] for x in vals], markersize=2)
	#scatter!(p, [x[1] for x in accurate_envelope], [x[2] for x in accurate_envelope])
	#display(p)
	return E
end

#=
function scatter_envelope!(p, E::Envelope)
	vals = [x[1] for x in E.A]
	@assert length(vals) > 0
	scatter!(p, [[x[i] for x in vals] for i in 1:length(vals[1])]..., markersize=2)
end
=#


function longitude_to_candidate(bt, longitude)
	d=Dict{Junction, Union{Piecewise,Linear}}()
	for J in bt.junctions
		inputs = []
		outputs = []
		for i in 0:J.left_len -1
			track = bt.forwardfan[(i,J)]
			n = round(Int,longitude[track[1]])
			for k in 1:n
				push!(inputs, i + k/(n+1))
			end
		end
		for i in 0:J.right_len-1
			track = bt.backwardfan[(i,J)]
			n = round(Int,longitude[track[1]])
			for k in 1:n
				push!(outputs, i + k/(n+1))
			end
		end
		@assert length(inputs)==length(outputs)
		d[J] = piecewise(inputs./ J.left_len, outputs./J.right_len)
	end
	return Candidate(bt,d)
end

function plotjs(A::Vector{Envelope}; maxabs=25)
	PlotlyJS.plot([_plotjs(E; maxabs=maxabs) for E in A])
end

function plotjs(E::Envelope; maxabs=25)
	plotjs(Envelope[E]; maxabs=maxabs)
end

function _plotjs(E::Envelope; maxabs=25, fill=false)
	pts = [x[1] for x in E.A if maximum(abs.(x[1]))<= maxabs]
	dim = length(pts[1])


	if dim==2
		if fill
			sort!(pts, by=x->x[1])
			all_pts=sort(vcat(pts,[(pts[i][1], pts[i+1][2]) for i in 1:length(pts)-1]), by=x->(x[1],-x[2]))
			PlotlyJS.scatter(x=[x[1] for x in all_pts],y=[x[2] for x in all_pts], mode="markers", fill="tozeroy")
			#construct the sequence
		else
			PlotlyJS.scatter(x=[x[1] for x in pts],y=[x[2] for x in pts], mode="markers")
		end
	elseif dim==3
		PlotlyJS.scatter(x=[x[1] for x in pts],y=[x[2] for x in pts], z=[x[3] for x in pts], mode="markers", type="scatter3d")
	else
		@assert false	
	end
end

function _plotjs(E1::Envelope{Lower}, E2::Envelope{Upper})
	pts1 = [x[1] for x in E1.A]
	pts2 = [x[1] for x in E2.A]

	@assert length(pts1[1])==2

	sort!(pts1, by=x->x[1])
	sort!(pts2, by=x->x[1])
	all_pts1=sort(vcat(pts1,[(pts1[i+1][1], pts1[i][2]) for i in 1:length(pts1)-1]), by=x->(x[1],-x[2]))
	all_pts2=sort(vcat(pts2,[(pts2[i][1], pts2[i+1][2]) for i in 1:length(pts2)-1]), by=x->(x[1],-x[2]))

	if all(all_pts1[1] .< all_pts2[1])
		pushfirst!(all_pts2, [all_pts1[1][1], all_pts2[1][2]])
	end

	if all(all_pts1[end] .< all_pts2[end])
		push!(all_pts1, [all_pts2[end][1], all_pts1[end][2]])
	end

	return [PlotlyJS.scatter(x=[x[1] for x in all_pts1],y=[x[2] for x in all_pts1], mode="markers", fill="tonexty", fillcolor="#00000000", marker=attr(color="#1f77b4")),
			PlotlyJS.scatter(x=[x[1] for x in all_pts2],y=[x[2] for x in all_pts2], mode="markers", fill="tonexty",marker=attr(color="#1f77b4")),
			]
end

function clear()
	deletetraces!(p,0:10)
end

function constraints(L::Longitude)
	ss = slopes(L)
	chi = -sum(L.weights)//2
	#@show chi

	npunctures = [gcd(a,b) for (a,b) in ss]
	#@show npunctures
	closed_chi = chi + npunctures[2] + npunctures[1]

	#=
	if closed_chi >= 0
		#@show L
		@show ss
		@show closed_chi
	end
	=#

	q,p = ss[1]
	s,r = ss[2]
	closed_chi = chi+npunctures[2]
	#=
		if q-closed_chi == 0 
			return (0,0)
		else
			return (p/(q - closed_chi),  r/s)
		end
	=#
	return [ (p/q,  r/(s+ (chi+npunctures[1]))),
			(p/(q+(chi+npunctures[2])), r/s)]


	#r/s, r/(s+chi)

	#(r+s)/s, (r+s)/(s+chi)
end

