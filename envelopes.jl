struct Envelope{S,T,D} #keep track of local maxes
    A::Vector{Tuple{Vector{T},D}}
	L::SpinLock
end

function Envelope{S,T,D}() where {S <: Comp, T, D}
	return Envelope{S,T,D}(Tuple{Vector{T},D}[], SpinLock())
end

function Envelope{S}(A::Vector{Tuple{Vector{T},D}}) where {S<: Comp, T, D}
	return Envelope{S,T,D}(A, SpinLock())
end

function Envelope()
    return Envelope{Upper,Float64,Cand{DiscreteHomeo}}()
end

function PEnvelope()
    return Envelope{Eq,Float64,Cand{DiscreteHomeo}}()
end

function strict_compare(x::Vector{T},y::Vector{S}) where {S,T}
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

#defining crevices
#order all points by x,y, and z coordinate.
#a crevice has the property that it is not in the interior of envelope, but if you move slightly down, left, or right, you enter the envelope. So it is a triple of points in the envelope, p(_x, p_y, p_z, such that 
#

#We'll start with the points [Inf, -Inf, -Inf], [-Inf, -Inf, Inf], [-Inf, Inf, -Inf].
#We'll maintain a tree-like structure, where each node is a triangle
#Whenever we insert a point, it restricts both the 
#

mutable struct Crevice{N,T} #can also be thought of as an octant in space
    faces::MVector{N,SVector{N,T}} #should be the points giving rise to this crevice
    pivot::Union{SVector{N,T},Nothing}
    children::Vector{Crevice{N}}
end
function Crevice(faces::AbstractVector{SVector{N,T}}) where {N,T}
    return Crevice{N,T}(MVector{N,SVector{N,T}}(faces), nothing, Crevice{N,T}[])
end

function contains(c::Crevice{N}, p::SVector{N,T}) where {N,T}
    all(c.faces[i][i] < p[i] for i in 1:N)
end

function is_valid_crevice(c::Crevice{N}) where {N}
    return all( (i==j || c.faces[i][i] < c.faces[j][i]+0.000001)
               for i in 1:N, j in 1:N
              )
end

function push!(c::Crevice{N}, pt::SVector{N,T}) where {N,T}
    if contains(c,pt)
        if c.pivot == nothing
            c.pivot=pt
            for i in 1:N
                childfaces = copy(c.faces)
                childfaces[i] = pt
                child = Crevice(childfaces)
                if is_valid_crevice(child)
                    push!(c.children, child)
                end
            end 
        else
            for child in c.children
                push!(child, pt)
            end
        end
    end 
end

function leaves(c::Crevice{N}) where {N}
    if c.pivot == nothing
        return [SVector{N}([c.faces[i][i] for i in 1:N])]
    else
        return Iterators.flatten(map(leaves, c.children))
    end
end

function all_leaves(c::Crevice{N}) where {N}
    return unique(collect(leaves(c)))
end

function crevices_general(e::Envelope{Upper,T}) where {T}
    N=length(e.A[1][1])
    c=Crevice([SVector{N,T}([(i==j ? -CLIP : CLIP) for j in 1:N]) for i in 1:N])
    for (x,_) in e.A
        push!(c, SVector{N,T}(x))
    end
    return [Vector{T}(x) for x in all_leaves(c)]
end

function crevices_general(e::Envelope{Lower,T}) where {T}
    N=length(e.A[1][1])
    c=Crevice([SVector{N,T}([(i==j ? -CLIP : CLIP) for j in 1:N]) for i in 1:N])
    for (x,_) in e.A
        push!(c, -SVector{N,T}(x))
    end
    return [Vector{T}(-x) for x in all_leaves(c)]
end


function inclosure(e::Envelope{S}, x::Vector{T}) where {S, T}
    return any(comp(S, x, y[1]) for y in e.A)
end

function push!(e::Envelope{S,T,D}, x::Tuple{Vector{T},D}) where {S,T,D}
	lock(e.L) do
		if !hasnan(x[1]) && !any(comp(S, x[1], y[1]) for y in e.A)
			filter!(y->!comp(S,y[1],x[1]), e.A)
			push!(e.A, x)
		end
	end
end

function push!(e::Envelope, e2::Envelope)
	for i in e2.A
		push!(e,i)
	end
end

function push!(e::Envelope{S,T,D}, _x::Tuple{Union{NTuple{N,R}, Vector{R}},D}) where {S, N, R <: Real,D,T}
    x = (T[_x[1]...], _x[2])
	lock(e.L) do
		if !hasnan(x[1]) && !any(comp(S, x[1], y[1]) for y in e.A)
			filter!(y->!comp(S,y[1],x[1]), e.A)
			push!(e.A, x)
		end
	end
end
