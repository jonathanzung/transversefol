module Envelopes
using Base.Threads
using StaticArrays
using ProgressMeter
export Upper, Lower, Eq, Envelope, crevices_general, basis_change
import Base: inv, getindex, setindex!, hash, push!, length, copy, show, rand

abstract type Comp
end

struct Upper <: Comp
end

struct Lower <: Comp
end

struct Eq <: Comp
end



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

length(e::Envelope) = length(e.A)

function hasnan(s)
	ret = any(map(isnan,s))
	if ret
		#println("rejecting $(s)")
	end
	return ret
end

function all_leq(x::Vector{T},y::Vector{S}) where {S,T}
	return all(x .<= y)
end

function all_lt(x::Vector{T},y::Vector{S}) where {S,T}
	return all(x .< y)
end
function strict_comp(::Type{Upper}, x::T, y::S) where {S,T}
	return x < y
end
function strict_comp(::Type{Lower}, x::T, y::S) where {S,T}
	return x > y
end

function comp(S::Type{Upper}, x, y)
	return all_leq(x, y)
end

function comp(S::Type{Lower}, x, y)
	return all_leq(y, x)
end

function comp(S::Type{Eq}, x, y)
	return x==y
end

function strict_comp(::Type{Upper}, x::Vector, y::Vector)
	return all_lt(x, y)
end

function strict_comp(::Type{Lower}, x::Vector, y::Vector)
	return all_lt(y, x)
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

function crevices_general(e::Envelope{Upper,T}, clip) where {T}
    N=length(e.A[1][1])
    c=Crevice([SVector{N,T}([(i==j ? -clip : clip) for j in 1:N]) for i in 1:N])
    @showprogress desc="computing crevices" for (x,_) in e.A
        push!(c, SVector{N,T}(x))
    end
    return [Vector{T}(x) for x in all_leaves(c)]
end

function crevices_general(e::Envelope{Lower,T}, clip) where {T}
    N=length(e.A[1][1])
    c=Crevice([SVector{N,T}([(i==j ? -clip : clip) for j in 1:N]) for i in 1:N])
    @showprogress desc="computing crevices" for (x,_) in e.A
        push!(c, -SVector{N,T}(x))
    end
    return [Vector{T}(-x) for x in all_leaves(c)]
end


function inclosure(e::Envelope{S}, x::Vector{T}) where {S, T}
    return any(comp(S, x, y[1]) for y in e.A)
end

function ininterior(e::Envelope{S}, x::Vector{T}) where {S, T}
    return any(strict_comp(S, x, y[1]) for y in e.A)
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

function rand!(e::Envelope)
    lock(e.L) do
        return rand(e.A)
    end
end


#Given lower_pts and upper_pts, compute a set of cubes whose union is the box closure

#=
function box_closure(lower_pts::Vector{SVector{N,T}}, upper_pts::Vector{SVector{N,T}}) where {N,T}
    if length(lower_pts) <= 1 || length(upper_pts)<=1
        return [(lp, up) for lp in lower_pts, up in upper_pts if strict_compare(lp,up)]
    end
    lower_pivot = rand(lower_pts)
    upper_pivot = rand(upper_pts)

    upper_splits_high = [filter(pt->pt[i] >= upper_pivot[i]),
    upper_splits_low = [filter(pt->pt[i] <= lower_pivot[i]) for i in 1:N]]

    lower_splits_high = [filter(pt->pt[i] >= upper_pivot[i]),
    lower_splits_low = [filter(pt->pt[i] <= lower_pivot[i]) for i in 1:N]]

    #Split into four batches:
    #upper_high, upper_low
    #lower_high, lower_low.

    #For points in upper

    return cat(
    [(lower_pivot,upper_pivot)],

    dims=1)


end

=#

# Apply a per-coordinate change of basis to the region (Elower, Eupper).
# transforms[i] = [a b; c d] with det=1; acts on slope x as f_i(x) = (c + d*x)/(a + b*x).
# The pole of f_i is -a/b (absent when b=0). Splits into up to 2^n pieces, one per sign
# pattern of each coordinate relative to its pole. Empty pieces are dropped.
function basis_change(Elower::Envelope{Lower,T,D}, Eupper::Envelope{Upper,T,D},
                      transforms::Vector{M}) where {T<:Rational, D, M<:AbstractMatrix}
    n = length(transforms)

    # Pole for each coordinate: Rational{Int} or nothing (b==0, affine, no split)
    poles = [transforms[i][1,2] == 0 ? nothing : -transforms[i][1,1]//transforms[i][1,2]
             for i in 1:n]

    # Möbius transformation for coordinate i
    f(i, x) = let a=transforms[i][1,1], b=transforms[i][1,2],
                  c=transforms[i][2,1], d=transforms[i][2,2]
        (c + d*x) / (a + b*x)
    end

    transform_s(s) = T[f(i, s[i]) for i in 1:n]

    # Choices per coordinate: nothing means no split, true/false = above/below pole
    choices = [poles[i] === nothing ? [nothing] : [true, false] for i in 1:n]

    result = Tuple{Envelope{Lower,T,D}, Envelope{Upper,T,D}}[]

    for pattern in Iterators.product(choices...)
        filter_fn(s) = all(
            poles[i] === nothing ||
            (pattern[i] ? s[i] > poles[i] : s[i] < poles[i])
            for i in 1:n)

        lower_pts = [(transform_s(s), c) for (s, c) in Elower.A if filter_fn(s)]
        upper_pts = [(transform_s(s), c) for (s, c) in Eupper.A if filter_fn(s)]

        if isempty(lower_pts) && isempty(upper_pts)
            continue
        end

        # When one side is empty, synthesize a boundary point from the pole images.
        # For coord i with a pole and pattern[i]=true:  f maps (pole,+∞) → (-∞, d/b),
        #   so the lower boundary is -∞ and the upper boundary is d/b.
        # For coord i with a pole and pattern[i]=false: f maps (-∞,pole) → (d/b,+∞),
        #   so the lower boundary is d/b and the upper boundary is +∞.
        # For affine coords (no pole): lower boundary is -∞, upper is +∞.
        if isempty(lower_pts)
            bound = T[if poles[i] === nothing
                          T(-1, 0)
                      elseif pattern[i]
                          T(-1, 0)
                      else
                          transforms[i][2,2] // transforms[i][1,2]
                      end for i in 1:n]
            lower_pts = [(bound, upper_pts[1][2])]
        elseif isempty(upper_pts)
            bound = T[if poles[i] === nothing
                          T(1, 0)
                      elseif pattern[i]
                          transforms[i][2,2] // transforms[i][1,2]
                      else
                          T(1, 0)
                      end for i in 1:n]
            upper_pts = [(bound, lower_pts[1][2])]
        end

        push!(result, (Envelope{Lower}(lower_pts), Envelope{Upper}(upper_pts)))
    end

    return result
end

function slice(v::Vector{T}, fillings::Vector{Tuple{Int,Int}}) where {T<:Real}
    return T[v[i] for i in 1:length(v) if fillings[i] == (0,0)]
end

#return a slice of an envelope with given filling slopes
function slice(E::Envelope{S,T,D}, fillings::Vector{Tuple{Int,Int}}) where {S, T<:Real, D}
    Eslice = Envelope{S,T,D}()
    for (s, c) in E.A
        if all(fillings[i] == (0,0) || strict_comp(S, fillings[i][2]//fillings[i][1],s[i]) for i in 1:length(fillings))
            push!(Eslice, (slice(s, fillings),c))
        end
    end
    return Eslice
end

end
