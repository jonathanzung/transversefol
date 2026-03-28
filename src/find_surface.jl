using LinearAlgebra
using JuMP, HiGHS, GLPK
using Base.Iterators 
using Random
using Base.Threads
#import Nemo
using OffsetArrays
using .MyLinearAlgebra
using LinearAlgebraX
using DataStructures
using ProgressMeter

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

# Like slopes(Longitude), but works for any integer weight vector (including negative entries).
# Returns a vector of Slope (one per cusp).
function multislope_vec(bt, weights::OffsetArray)
    [sum(weights[i]*bt.weights[(i,j)] for (i,j) in edgelist) for edgelist in bt.alledges]
end

function fan_matrix(fans)
    relations = [((x[1] for x in f1), (x[1] for x in f2)) for (f1,f2) in fans]
    l = collect(Iterators.flatten(Iterators.flatten(relations)))
    n = maximum(l) + 1
    M = zeros(Int, length(relations), n)
    for (i, (l1, l2)) in enumerate(relations)
        for j in l1
            M[i, j+1] += 1
        end
        for j in l2
            M[i, j+1] -= 1
        end
    end
    return M
end

# Boundary operator d3: C3 -> C2 in the Regina convention.
# tet_faces[i] = [(tri_idx, sign), ...] for the 4 faces of tetrahedron i.
function d3_matrix(tet_faces, ntri)
    ntet = length(tet_faces)
    D = zeros(Int, ntri, ntet)
    for (j, faces) in enumerate(tet_faces)
        for (tri_idx, s) in faces
            D[tri_idx + 1, j] += s
        end
    end
    return D
end

# Returns a Z-basis for H_2(M, dM; Z) as OffsetArrays (weight vectors indexed from 0).
# H_2(M, dM) = ker(fan_matrix) / im(d3), where im(d3) <= ker(fan_matrix) by the switch conditions.
# face_coorientations is used to convert d3 from Regina convention to veering convention.
function compute_H2_rel_boundary(fans, tet_faces, face_coorientations)
    M = fan_matrix(fans)
    ntri = size(M, 2)

    # d3 in veering convention: multiply each row by face_coorientations
    D = d3_matrix(tet_faces, ntri)
    for i in 1:ntri
        D[i, :] .*= face_coorientations[i - 1]
    end

    # In a veering triangulation each tetrahedron has exactly 2 top and 2 bottom faces.
    for (j, faces) in enumerate(tet_faces)
        veering_signs = [face_coorientations[tri_idx] * s for (tri_idx, s) in faces]
        npos = count(>(0), veering_signs)
        nneg = count(<(0), veering_signs)
        @assert npos == 2 && nneg == 2 "Tetrahedron $j has $npos top and $nneg bottom faces (expected 2 each)"
    end

    @assert all(M * D .== 0) "im(d3) is not contained in ker(fan_matrix)"

    # Find integer Z-basis for ker(M) using Smith normal form of M.
    # M = L_M * S_M * U_M  =>  ker(M) = U_M^{-1} * { w : S_M*w = 0 }
    _, S_M, _, _, Uinv_M = MyLinearAlgebra.LDU(M)
    # Column j of Uinv_M is in ker(M) iff S_M[:,j] == 0.
    # S_M is (m×n); for j > m all columns are zero; for j ≤ m check the diagonal.
    m_M = min(size(S_M)...)
    ker_mask = [i > m_M || iszero(S_M[i,i]) for i in 1:size(Uinv_M, 2)]
    K = Uinv_M[:, ker_mask]  # integer Z-basis for ker(M), shape ntri × dim(ker)

    @assert all(M * K .== 0) "K is not in ker(M)"

    # Express im(D) within ker(M): solve K * A = D for integer A.
    # Use rational arithmetic to get the exact solution.
    K_rat = Rational{Int}.(K)
    A_rat = K_rat \ Rational{Int}.(D)
    A_int = round.(Int, A_rat)
    @assert K * A_int == D "d3 image does not lie exactly in ker(fan_matrix) basis"

    # H_2(M, dM) = cokernel of A: Z^ntet -> Z^dim(ker) / im(A).
    # Smith decomposition: A = L * S * U. Generators of Z^dim(ker) / im(S) lift to columns of L.
    L, S, _, _, _ = MyLinearAlgebra.LDU(A_int)
    r = min(size(A_int)...)
    torsion = [S[j,j] for j in 1:r if abs(S[j,j]) > 1]
    if !isempty(torsion)
        @warn "H_2(M, dM) has torsion: $torsion"
    end
    free_cols = [j for j in 1:size(A_int, 1) if j > r || S[j,j] == 0]
    return [OffsetArrays.Origin(0)(K * L[:, j]) for j in free_cols]
end

# Returns all simple directed cycles in the dual graph as vectors of edge indices
# (1-indexed into top_bot_pairs). Each cycle is reported exactly once, starting
# from its minimum vertex, so no duplicates arise from different starting points.
function all_simple_cycles(top_bot_pairs)
    V = maximum(Iterators.flatten(top_bot_pairs)) + 1
    # adj[v+1]: outgoing edges from vertex v (0-indexed), as (w+1, edge_idx) pairs
    adj = [Tuple{Int,Int}[] for _ in 1:V]
    for (i, (u, v)) in enumerate(top_bot_pairs)
        push!(adj[u+1], (v+1, i))
    end

    cycles = Vector{Vector{Int}}()

    function dfs(s, v, visited, path_edges)
        for (w, e) in adj[v]
            if w == s
                push!(cycles, push!(copy(path_edges), e))
            elseif w > s && w ∉ visited
                push!(visited, w)
                push!(path_edges, e)
                dfs(s, w, visited, path_edges)
                pop!(path_edges)
                delete!(visited, w)
            end
        end
    end

    for s in 1:V
        dfs(s, s, Set{Int}([s]), Int[])
    end

    return cycles
end

# Find all integer coordinate vectors a ∈ ℤᵇ¹ lying in the cone
#   { a : ⟨γ, Σₖ aₖ gₖ⟩ ≥ 0  for all dual-graph cycles γ }
# where the pairing of cycle γ with weight vector w is Σₑ∈γ w[e-1]
# (e is 1-indexed into top_bot_pairs; h2_gens are OffsetArrays indexed from 0).
# Returns only vectors with w·a < max_weight, where w[k] = sum(h2_gens[k]) is the pairing
# with the cycle that intersects each face exactly once.
function cone_integer_points(cycles, h2_gens; max_weight=20)
    b1 = length(h2_gens)
    b1 == 0 && return [Int[]]

    # P[j,k] = pairing of cycle j with generator k
    P = [sum(h2_gens[k][e-1] for e in cycle; init=0) for cycle in cycles, k in 1:b1]

    # w[k] = total weight of generator k = pairing with the cycle that intersects each face once
    w = [sum(h2_gens[k]) for k in 1:b1]

    # Use LP to compute tight per-dimension bounds: for each coordinate k, find the max/min value
    # of a_k over the polytope { w·a <= max_weight-1, P*a >= 0 }.
    # The LP is bounded because any unbounded ray d would require P*d >= 0 and
    # w·d <= 0 simultaneously. The cycle constraints collectively prevent this:
    # a direction that increases a_k while decreasing others (keeping weight bounded)
    # must eventually pair negatively with some cycle.
    lo = zeros(Int, b1)
    hi = zeros(Int, b1)
    for k in 1:b1
        for (sense, store) in ((MAX_SENSE, hi), (MIN_SENSE, lo))
            model = Model(HiGHS.Optimizer)
            set_silent(model)
            @variable(model, a[1:b1])
            @constraint(model, dot(w, a) <= max_weight - 1)
            @constraint(model, P * a .>= 0)
            @objective(model, sense, a[k])
            optimize!(model)
            @assert termination_status(model) == OPTIMAL "LP for cutoff estimation is unbounded or infeasible (coordinate $k, sense $sense)"
            store[k] = (sense == MAX_SENSE ? floor : ceil)(Int, value(a[k]))
        end
    end

    results = Vector{Vector{Int}}()
    for a in Iterators.product((lo[k]:hi[k] for k in 1:b1)...)
        av = collect(Int, a)
        if dot(w, av) < max_weight && all(>=(0), P * av)
            push!(results, av)
        end
    end
    return results
end

#=
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
=#

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

#compute a basis for H_1(M), and also try to include [1,1,....,1] in the basis
#Note that [1,1,...,1] is always a cycle, a union of branch loops. It might be homologically trivial though, for example in the case where the pseudo-Anosov flow admits no transverse surface.
function compute_homology(fans, top_bot_pairs)

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

function find_longitudes_hom2(fans, top_bot_pairs, tet_faces, face_coorientations;
                               primitive_only=true, max_weight=50)
    h2_gens = compute_H2_rel_boundary(fans, tet_faces, face_coorientations)
	#@show length(h2_gens)
    b1 = length(h2_gens)
    b1 == 0 && return OffsetArray[]

    cycles = all_simple_cycles(top_bot_pairs)
    @assert length(cycles) > 0 "no simple cycles found in dual graph"
    candidates = cone_integer_points(cycles, h2_gens; max_weight=max_weight)

    ntri = length(first(h2_gens))

    # d3 matrix in veering convention (same computation as compute_H2_rel_boundary)
    D = d3_matrix(tet_faces, ntri)
    for i in 1:ntri
        D[i, :] .*= face_coorientations[i-1]
    end
    ntet = size(D, 2)

    M = fan_matrix(fans)
    @assert all(M * D .== 0) "im(D) does not satisfy fan balance: M*D != 0"
    for k in 1:b1
        @assert all(M * collect(h2_gens[k]) .== 0) "h2_gens[$k] does not satisfy fan balance"
    end

    results = Vector{OffsetVector{Int}}(undef, length(candidates))

    @showprogress desc="Computing longitudes" @threads for idx in eachindex(candidates)
        a = candidates[idx]
        all(iszero, a) && continue

        # Base weight vector for this H₂ class
        h = sum(a[k] * collect(h2_gens[k]) for k in 1:b1)
        @assert all(M * h .== 0) "h does not satisfy fan balance for a=$a"

        # Find t ∈ ℤ^ntet minimising sum(x) subject to x = h + D*t ≥ 0
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        @variable(model, t[1:ntet], Int)
        x_expr = D * t .+ h
        @constraint(model, x_expr .>= 0)
        @constraint(model, sum(x_expr) >= 1)
        @objective(model, Min, sum(x_expr))
        optimize!(model)

        @assert is_solved_and_feasible(model) "LP should always be feasible: every H₂ class in the positive cone has a non-negative representative"
        t_val = round.(Int, value.(t))
        x_val = OffsetArrays.Origin(0)(D * t_val .+ h)
        @assert all(M * collect(x_val) .== 0) "solution x_val does not satisfy fan balance"
        @assert eltype(x_val) == Int "x_val has unexpected element type $(eltype(x_val))"
        if !primitive_only || gcd(a...) == 1
            @assert all(sum(x_val[e-1] for e in cycle; init=0) >= 0 for cycle in cycles) "longitude does not pair non-negatively with all cycles"
            results[idx] = x_val
        end
    end

    return sort!([results[i] for i in eachindex(results) if isassigned(results, i)], by=x->sum(x))
end


#=
function find_longitudes_hom(fans, top_bot_pairs; primitive_only=true) #find longitudes by homology class
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
                if !primitive_only || is_primitive(O)
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
=#

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
