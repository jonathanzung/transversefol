using LinearAlgebra
using JuMP, HiGHS, GLPK
using Base.Iterators 
using Random
using Base.Threads
#import Nemo
using OffsetArrays
include("MyLinearAlgebra.jl")
using .MyLinearAlgebra
using LinearAlgebraX
using DataStructures

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


#It would be better to first compute the 

function find_longitude(fans)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(Iterators.flatten(Iterators.flatten(relations)))
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
	#model = Model(Cbc.Optimizer)
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

function connected_components(l, fans)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

    faces = Iterators.flatten([(i,j) for j in 1:l[i]] for i in eachindex(l)) |>  collect
    A=DisjointSets(faces)

    for (r1,r2) in relations
        stack1 = Iterators.flatten([[(i,j) for j in 1:l[i]] for i in r1]) |> collect
        stack2 = Iterators.flatten([[(i,j) for j in 1:l[i]] for i in r2]) |> collect
        @assert length(stack1) == length(stack2)
        for (x1,x2) in zip(stack1,stack2)
            union!(A,x1,x2)
        end
    end
    return length(unique([find_root(A,x) for x in faces]))
end

function is_fiber(l, top_bot_pairs)#check if a longitude is a fiber surface. We do this by checking if every cycle in the dual graph intersects the fiber.
    uncut_edges = [top_bot_pairs[i] for i in 1:length(top_bot_pairs) if l[i-1]==0] 

    #now check if this graph has any cycles
    
    V=maximum(Iterators.flatten(top_bot_pairs))+1

    @label here
    vert_out_degrees=OffsetArrays.Origin(0)([0 for i in 1:V])
    for (i,j) in uncut_edges
        vert_out_degrees[i]+=1 
    end
    N=length(uncut_edges)
    for i in eachindex(vert_out_degrees)
        if vert_out_degrees[i]==0
            filter!(x->x[2] != i, uncut_edges)
        end
    end

    if length(uncut_edges) < N
        @goto here
    end

    return length(uncut_edges) == 0
end
function compute_homology(fans, top_bot_pairs)
    #compute a basis for H_1, and also try to include [1,1,....,1] in the basis
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

    #V,E,F are the number of vertices, edges, faces in the dual graph
	l = collect(Iterators.flatten(Iterators.flatten(relations)))
	sort!(l)
	E = maximum(l)+1 #equals number of faces in the veering triangulation
    F=length(relations) #equals number of edges in the veering triangulation
	@assert length(l) == (maximum(l)+1)*3

    l2 = unique(collect(Iterators.flatten(top_bot_pairs)))
    V = maximum(l2)+1

    @assert V-E+F==0
    

    M=zeros(Int, V+E+F, V+E+F)

    for (f,r) in enumerate(relations)
        for i in r[1]
            M[V+i+1,V+E+f]+=1
        end
        for i in r[2]
            M[V+i+1,V+E+f]-=1
        end
    end

    for (e,r) in enumerate(top_bot_pairs)
        M[r[1]+1, V+e]+=1
        M[r[2]+1, V+e]-=1
    end

    diag, gens = MyLinearAlgebra.homology(M)
    allgens = [gens[:,i] for i in 1:size(gens,2) if diag[i]==0] #extract the generators corresponding to nontorsion
    #vertgens = filter(x->all(x[V+1:E].==0) && all(x[V+E+1:V+E+F].==0), allgens) #extract the generators of H_0
    edgegens = filter(x->all(x[1:V].==0) && all(x[V+E+1:V+E+F].==0), allgens) #extract the generators of H_1
    #facegens = filter(x->all(x[1:V+E].==0), allgens) #extract the generators of H_2
    trimmed = [x[V+1:V+E] for x in edgegens] #the generators of H_1, as a vector of length E.

    extra = [1 for i in 1:E] #we want this element in our homology basis. So we will add it first, and then remove elements which are dependent on previous ones

    pushfirst!(trimmed, extra)

    d_image = [M[V+1:V+E, i] for i in V+E+1:V+E+F]

    lastrank = rankx(stack(d_image)) #a basis for the image of d

    @show length(trimmed)
    @show lastrank

    for i in 1:length(trimmed)
        @show trimmed

        nextrank = rankx(stack(cat(d_image,trimmed[1:i],dims=1)))
        @show nextrank
        if nextrank == lastrank
            if i==1
                println("Warning, all 1's is trivial in homology")
            end
            deleteat!(trimmed, i)
            break
        end
        lastrank = nextrank
    end
    @show trimmed
    return trimmed 
end

function find_longitudes_hom(fans, top_bot_pairs) #find longitudes by homology class
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(Iterators.flatten(Iterators.flatten(relations)))
	sort!(l)
	n = maximum(l)+1
	@assert length(l) == (maximum(l)+1)*3
	#println(n)

    function indicator(hom_class)
        A=zeros(Int, n)
        for i in hom_class
            A[i+1]+=1
        end
        return A
    end

    hom_classes = compute_homology(fans, top_bot_pairs)
    @assert length(hom_classes) >= 1

    H = transpose(stack(hom_classes))
    nH = size(H,1)
    H=collect(H)
    @show collect(H)
    #@show typeof(collect(H))


	M=zeros(length(relations),n)
	for (i,(l1,l2)) in enumerate(relations)
		for j in l1
			M[i,j+1] += 1
		end
		for j in l2
			M[i,j+1] -= 1
		end
	end

	ch=Channel(3*Threads.nthreads())

	function search_interval(a)
        println("Searching $(a) from $(Threads.threadid())")
		model = Model(HiGHS.Optimizer)
		set_silent(model)
		@variable(model, x[1:n], Int)
		@constraint(model, M * x .== 0)
		@constraint(model, x .>= 0)
		#@constraint(model, sum(x) == a)

		#@objective(model, Min, sum(x))


        @constraint(model, con, H*x == [0 for i in 1:nH])

        for hom_class in CartesianIndices(tuple([-3*a*sum(abs.(hom_classes[i])):3*a*sum(abs.(hom_classes[i])) for i in 2:nH]...))
            if is_valid(model, con)
                delete(model, con)
                unregister(model, :con)
            end
            @constraint(model, con, H*x == Int[a, Tuple(hom_class)...])
            @objective(model, Min, 0)
            optimize!(model)
            if is_solved_and_feasible(model)
                val=round.(Int,value.(x))
                O=OffsetArrays.Origin(0)(val)
                if is_primitive(O)
                    try
                        put!(ch, O)
                    catch e
                        if e isa InvalidStateException
                            #println("terminating")
                            break
                        else
                            rethrow(e)
                        end
                    end
                end
            end
            if !isopen(ch)
                #println("terminating")
                break
            end
		end
	end

    @async begin
        @threads :greedy for a in takewhile(x -> isopen(ch), 1:300) #Iterators.countfrom(0,1))
            search_interval(a)
        end
        close(ch)
    end

	#todo: clean up tasks when channel is closed
	return ch

end

#=
function find_longitudes_random(fans)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(Iterators.flatten(Iterators.flatten(relations)))
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

	ch=Channel(10*Threads.nthreads())

	function search_interval(a)
		println("Searching $((a))")
		model = Model(HiGHS.Optimizer)
		set_silent(model)
		@variable(model, x[1:n], Int)
		@constraint(model, M * x .== 0)
		@constraint(model, x .>= 0)
		@constraint(model, sum(x) == a)

		#@objective(model, Min, sum(x))


         @constraint(model, con, sum(x) >= 0)
		for i in 1:50
			perturbation = randn(n)
            if is_valid(model, con)
                delete(model, con)
                unregister(model, :con)
            end
            @constraint(model, con, sum(x.*perturbation) >= 0)
			@objective(model, Min, sum(x .* perturbation)) #problem: this prefers ones with larger L^2 norm
			optimize!(model)
			if is_solved_and_feasible(model)
				val=round.(Int,value.(x))
				O=OffsetArrays.Origin(0)(val)
				if is_primitive(O)
					put!(ch, O)
				end
			else
				break
			end
		end
	end

	@async @threads :greedy for a in Iterators.countfrom(0,1)
		search_interval(a)
	end

	#todo: clean up tasks when channel is closed
	return ch

end

function find_longitudes_iterative(fans, maxchi)
	relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]

	l = collect(Iterators.flatten(Iterators.flatten(relations)))
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

	perturbation = 1 .+ 0.1 .* rand(Xoshiro(123), n)
	ch=Channel(10*Threads.nthreads())

	function search_interval(a,b)
		println("Searching $((a,b))")
		model = Model(HiGHS.Optimizer)
		set_silent(model)
		@variable(model, x[1:n], Int)
		@constraint(model, M * x .== 0)
		@constraint(model, x .>= 0)
		@constraint(model, sum(x) >= 1)
		@constraint(model, sum(x .* perturbation) >= a)

		@objective(model, Min, sum(x .* perturbation))
		#@objective(model, Min, sum(x))


		while true
			optimize!(model)
			val=round.(Int,value.(x))
			objval=sum(value.(x) .* perturbation)
			O=OffsetArrays.Origin(0)(val)
			if is_primitive(O)
				put!(ch, O)
			end
			if objval > b
				break
			end
			@constraint(model, sum(x .* perturbation) >= objval + 0.001)
			#println(collect((i-1)=>j for (i,j) in enumerate(value.(x)) if j!=0))
			#return Dict((i-1)=>j for (i,j) in enumerate(value.(x)))
		end
	end

	step=2
	@async @threads :greedy for a in Iterators.countfrom(0,step)
		search_interval(a,a+step)
		if a + step > maxchi
			break
		end
	end

	#todo: clean up tasks when channel is closed
	return ch

end
=#

#=
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
=#
