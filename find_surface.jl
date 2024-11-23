using LinearAlgebra
using JuMP, HiGHS
using Base.Iterators 
import Nemo
using OffsetArrays

#=
relations10_145 = 
[
 ((3,5,9,8,7),(6,1,0,2,4)),
 ((4,10),(3,11)),
 ((11,7),(10,6)),
 ((1,7,10,3,2),(0,0)),
 ((1,2),(8,5)),
 ((5,4,11,6,8),(9,9))
]
relations10_139 = 
[
 ((4,7,6,0),(1,3,1)),
 ((0,2,5,4),(9,8,9)),
 ((6,9,7),(5,1,2)),
 ((2,6),(3,4,8)),
 ((8,0,3),(7,5))]
 =#

function is_primitive(A)
	return gcd(A...)==1
end

function find_longitude(fans)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(flatten(flatten(relations)))
	sort!(l)
	n = maximum(l)+1
	@assert length(l) == (maximum(l)+1)*3
	#println(n)



	M=zeros(length(relations),n)
	for (i,(l1,l2)) in enumerate(relations)
		for j in l1
			M[i,j+1] += 1
		end
		for j in l2
			M[i,j+1] -= 1
		end
	end

	model = Model(HiGHS.Optimizer)
	set_silent(model)
	@variable(model, x[1:n], Int)
	@constraint(model, M * x .== 0)
	@constraint(model, x .>= 0)
	@constraint(model, sum(x) >= 1)
	@objective(model, Min, sum(x))
	optimize!(model)
	#println(collect((i-1)=>j for (i,j) in enumerate(value.(x)) if j!=0))
	return OffsetArrays.Origin(0)(value.(x))
	#return Dict((i-1)=>j for (i,j) in enumerate(value.(x)))
end

function find_longitudes_iterative(fans, nsolutions)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(flatten(flatten(relations)))
	sort!(l)
	n = maximum(l)+1
	@assert length(l) == (maximum(l)+1)*3
	#println(n)



	M=zeros(length(relations),n)
	for (i,(l1,l2)) in enumerate(relations)
		for j in l1
			M[i,j+1] += 1
		end
		for j in l2
			M[i,j+1] -= 1
		end
	end

	model = Model(HiGHS.Optimizer)
	set_silent(model)
	@variable(model, x[1:n], Int)
	@constraint(model, M * x .== 0)
	@constraint(model, x .>= 0)
	@constraint(model, sum(x) >= 1)
	perturbation = [1+0.01*rand() for i in 1:n]
	@objective(model, Min, sum(x .* perturbation))
	#@objective(model, Min, sum(x))

	ret = []
	while length(ret) < nsolutions
		optimize!(model)
		val=round.(Int,value.(x))
		objval=sum(value.(x) .* perturbation)
		O=OffsetArrays.Origin(0)(val)
		if is_primitive(O)
			push!(ret, O)
		end
		@constraint(model, sum(x .* perturbation) >= objval + 0.00001)
		#println(collect((i-1)=>j for (i,j) in enumerate(value.(x)) if j!=0))
		#return Dict((i-1)=>j for (i,j) in enumerate(value.(x)))
	end
	return ret
end

function find_longitudes(fans)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(flatten(flatten(relations)))
	sort!(l)
	n = maximum(l)+1
	@assert length(l) == (maximum(l)+1)*3
	
	M=zeros(length(relations),n)
	for (i,(l1,l2)) in enumerate(relations)
		for j in l1
			M[i,j+1] += 1
		end
		for j in l2
			M[i,j+1] -= 1
		end
	end
	S = Nemo.ZZMatrixSpace(length(relations),n)
	(d, bmat) = nullspace(S(M))
	#@show bmat, d
	ret = []
	for ind in CartesianIndices(tuple((-11:11 for i in 1:d)...))
		tmp = OffsetArrays.Origin(0)(Vector{Int}(bmat*[Tuple(ind)...]))
		if all(tmp .>= 0)
			push!(ret, tmp)
		end
	end
	return ret
end
