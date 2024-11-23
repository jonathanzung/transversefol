const T=Float64
using OffsetArrays
import Base: inv, getindex, setindex!, hash, push!, length
using DataStructures
using Plots
using Profile
using ProfileView
using Base.Threads
using PlotlyJS

plotlyjs()

import Plots: plot


abstract type Homeo end

#really, the identity function
struct Linear <: Homeo
end

#represents a piecewise linear function from [0,1] to [0,1]
struct Piecewise <: Homeo
	left::Union{Piecewise, Linear}
	right::Union{Piecewise, Linear}
	x::T
	fx::T #x maps to fx
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

struct Junction
	index::Int
	inv::Bool
	left_len::Int
	right_len::Int
end

function Base.hash(J::Junction, h::UInt)
	return hash(J.index, hash(J.inv, h))
end

const Track = Tuple{Int,Int} #edge index and vertex number

struct State
	x::T#height in the edge
	e::Track#which track we're on
end

function inv(e::Junction)
	return Junction(e.index,!e.inv,e.right_len,e.left_len)
end

struct BoundaryTriangulation
	forward::Dict{Track, Tuple{Int,Junction}}
	backward::Dict{Track, Tuple{Int,Junction}}
	forwardfan::Dict{Tuple{Int,Junction}, Track}
	backwardfan::Dict{Tuple{Int,Junction}, Track}
	weights::DefaultDict{Track, Vector{T}}
	junctions::Vector{Junction}
	firstrungs::Vector{Track}#one rung in each cusp
	#fans is a list of pairs of fans
	
	function BoundaryTriangulation(fans, face_coorientations, firstrungs)
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

		weights=DefaultDict{Track,Vector{T}}(T[0,0])
		return new(forward,backward,forwardfan,backwardfan,weights,junctions,firstrungs)
	end
end

struct Candidate
	bt::BoundaryTriangulation
	d::Dict{Junction, Union{Piecewise,Linear}}
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

	return State(fx-Int(floor(fx)), bt.backwardfan[(j,J)])
end

function trace_backwards(s::State, bt::BoundaryTriangulation, c::Candidate)
	i,J = bt.backward[s.e]
	fx = c(inv(J), i+s.x)
	j=Int(floor(fx))
	return State(fx-Int(floor(fx)), bt.forwardfan[(j,J)])
end




function annealing(f, initial, jiggle, betastart, betafinish, nsteps; verbose=false)
	#linear annealing on range betastart, betafinish
	current = initial
	currval = f(current)
	reject_count = 0
	accept_count = 0

	vals=[]
	push!(vals,currval)
	for (i,currbeta) in zip(1:nsteps, range(betastart,betafinish, nsteps))
		jig = jiggle(current)
		newval = f(jig)
		prob = exp(currbeta * newval)/( exp(currbeta * currval)  + exp(currbeta*newval))
		#@show prob
		if rand() < exp(currbeta * newval)/( exp(currbeta * currval)  + exp(currbeta*newval))
			current = jig
			currval = newval
			accept_count += 1
		else
			reject_count += 1
		end
		push!(vals, currval)

		if verbose && i%10000 == 0
			@show (all_slopes(current))
			@show (accept_count, reject_count) 
			p=plot(1:length(vals), vals)
			display(p)
		end

	end
	@show (accept_count, reject_count)
	return (vals,current)
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

function all_slopes(c::Candidate; time=200)
	return [slope(c; s=State(0.023423,rung), time=time) for rung in c.bt.firstrungs]
end

function slope(c::Candidate; time=200, s=State(0.312423, (1,0)))
	weight=T[0,0]
	weight .+= c.bt.weights[s.e]
	while abs(weight[1]) + abs(weight[2]) < time
		#println(s)	
		s=trace_forwards(s, c.bt, c)
		weight .+= c.bt.weights[s.e]
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
	return weight[2]/weight[1]
end

struct Envelope #keep track of local maxes
	A::Vector{Tuple{Vector{T},Candidate}}
	comp::Function
end

function Envelope()
	return Envelope(Tuple{Vector{T},Candidate}[], strict_compare)
end

function PEnvelope()
	return Envelope(Tuple{Vector{T},Candidate}[], (x,y)->false)
end

function strict_compare(x::Vector{T},y::Vector{T})
	return all(x .<= y)
end

function push!(e::Envelope,  x::Tuple{Vector{T},Candidate})
	if !any(strict_compare(x[1], y[1]) for y in e.A)
		filter!(y->!e.comp(y[1],x[1]), e.A)
		push!(e.A, x)
	end
end

function push!(e::Envelope, x::Candidate)
	push!(e, (all_slopes(x,time=5000), x))
end

function relu(x)
	if x > 0
		x
	else
		0
	end
end

function objective(slope, old_slope, target)
	return sum(relu.(slope.-old_slope) .- 10 * relu.(slope .- target) .- 10 * relu.(old_slope .- slope))
end

function try_improve(E::Envelope; nsubdivide=0, iters=50000, time=1000, target = [1000,1000])
	accurate_E=Envelope()
	l=SpinLock()
	@threads for i in 1:length(E.A)
		_,oldcand = E.A[i]
		old_v = all_slopes(oldcand,time=30000)

		for i in 1:nsubdivide
			oldcand = subdivide(oldcand)
		end
		lock(l) do
			push!(accurate_E, (old_v, oldcand))
		end
		if objective(old_v, old_v, target) < -0.001
			continue
		end

		_,newcand = annealing(c->objective(all_slopes(c,time=time), old_v, target), oldcand, c->jiggle(c,0.003), 400, 400, iters; verbose=false)
		new_v = all_slopes(newcand,time=30000)
		@show old_v, new_v
		lock(l) do
			push!(accurate_E, (new_v, newcand))
		end
	end
	return accurate_E
end

length(e::Envelope) = length(e.A)

function random_trials(bt)
	N=200000
	E=Envelope()
	l=SpinLock()
	trials=[T[] for i in 1:N]


	inrange(x) = all(-5 < i < 5 for i in x)

	for i in 1:N
		c=random_candidate(bt,1)
		trials[i] = all_slopes(c)
		lock(l) do
			push!(E, (trials[i],c))
		end
	end

	@show [x[1] for x in E.A]
	vals = collect(filter(inrange, trials))

	accurate_E = try_improve(E)
	accurate_envelope = sort(collect(filter(inrange, [v for (v,c) in accurate_E.A])))

	p=scatter([x[1] for x in vals], [x[2] for x in vals], markersize=2)
	scatter!(p, [x[1] for x in accurate_envelope], [x[2] for x in accurate_envelope])
	display(p)
end

function scatter_envelope!(p, E::Envelope)
	vals = [x[1] for x in E.A]
	@assert length(vals) > 0
	scatter!(p, [[x[i] for x in vals] for i in 1:length(vals[1])]..., markersize=2)
end


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

function plotjs(A::Vector{Envelope}; maxabs=5)
	PlotlyJS.plot([_plotjs(E; maxabs=maxabs) for E in A])
end

function plotjs(E::Envelope; maxabs=5)
	plotjs(Envelope[E]; maxabs=maxabs)
end

function _plotjs(E::Envelope; maxabs=5)
	pts = [x[1] for x in E.A if maximum(abs.(x[1]))<= maxabs]
	dim = length(pts[1])

	if dim==2
		PlotlyJS.scatter(x=[x[1] for x in pts],y=[x[2] for x in pts], mode="markers")
	elseif dim==3
		PlotlyJS.scatter(x=[x[1] for x in pts],y=[x[2] for x in pts], z=[x[3] for x in pts], mode="markers", type="scatter3d")
	else
		@assert false	
	end
end


include("find_surface.jl")
include("batch/manifest.txt")
isosig = "siddhi2"
isosig = "eLMkbcddddedde_2100"
#isosig = "gvLQQcdeffeffffaafa_201102"
begin
	println(isosig)
	include("batch/$(isosig).txt")
	bt=BoundaryTriangulation(fans, face_coorientations,firstrungs)
	for (x,y) in weights
		bt.weights[x] = y
	end
	#=
	global longitude
	if longitude == nothing
		longitude=find_longitude(fans)#weights of the different faces
	end
	=#
	longitudes = find_longitudes_iterative(fans,1000)
	#push!(longitudes, longitude)
	#sort!(longitudes, by=l->all_slopes(longitude_to_candidate(bt,l),time=30000)[1])
	function valid_slope(s)
		return maximum(abs.(s)) <= 5 && s[1] <= 0
	end
	filter!(l-> valid_slope(all_slopes(longitude_to_candidate(bt,l),time=2000)), longitudes)

	#c=random_candidate(bt,3)
	#=
	for i in 0:5
		for j in 0:2
			println(slope(c; time=10000, s = State(0.5, (i,j))) |> Float64)
		end
	end
	=#
	#
	
	#println(all_slopes(c, time=10000))
	
	dummy_candidate=random_candidate(bt,0)
	L6a2E=Envelope()
	for pt in [(-2,1/2),   (-1, 1/3) ,   (-1/2, 1/6),    (-1/3, 1/9),    (-1/6, 1/18)]
		push!(L6a2E, ([pt...], dummy_candidate))
	end

	if true
		E=Envelope()
		E2=PEnvelope()
		E3=PEnvelope()
		for l in longitudes
			local c
			xi = sum(l)
			c=longitude_to_candidate(bt,l)
			x=all_slopes(c, time=5000)
			push!(E, (x,c))
			#push!(E2, ([1/(1/x[1]-xi/2+1), x[2]], c))
			push!(E3, ([x[1], 1/(1/x[2]-xi/2+1)], c))
			println(x)
		end
		accurate_E = try_improve(E; nsubdivide=2, iters=100000, time=1000, target=[-1.5,0.5])

		p=plotjs([E,L6a2E,E3,accurate_E])
		display(p)
		#scatter_envelope!(p, accurate_E)
		#scatter!(p, [-2,-1,-1/2,0],[1/2,1/3,1/6,0])
		#scatter!(p,[-2,-1,-1/2],[1/3,1/6,0])
	end

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

end

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


