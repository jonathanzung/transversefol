module Envelopes
using Base.Threads
using StaticArrays
using ProgressMeter
export Upper, Lower, Eq, Envelope, crevices_general, basis_change, BasisChange
import Base: inv, getindex, setindex!, hash, push!, length, copy, show, rand, clamp, *

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
    return all( (i==j || c.faces[i][i] <= c.faces[j][i])
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

function accumulate_leaves!(c::Crevice{N},accum::Vector{SVector{N}}) where {N}
    if c.pivot == nothing
        push!(accum, SVector{N}([c.faces[i][i] for i in 1:N]))
    else
        for child in c.children
            accumulate_leaves!(child, accum)
        end
    end
    return accum
end

function all_leaves(c::Crevice{N}) where {N}
    return unique(accumulate_leaves!(c,SVector{N}[]))
end

function crevices_general(e::Envelope{Upper,T}, clip; N=length(e.A[1][1])) where {T}
    c=Crevice([SVector{N,T}([(i==j ? -clip : clip) for j in 1:N]) for i in 1:N])
    @showprogress desc="computing crevices" for (x,_) in e.A
        push!(c, SVector{N,T}(x))
    end
    return [Vector{T}(x) for x in all_leaves(c)]
end

function crevices_general(e::Envelope{Lower,T}, clip; N=length(e.A[1][1])) where {T}
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

# Möbius transform f(x) = (c+dx)/(a+bx) for matrix [a b; c d] with det=1.
# above_pole: true  (x ≥ pole, lo=pole) → pole maps to -∞
#             false (x ≤ pole, hi=pole) → pole maps to +∞
# For ±∞ input: affine (b=0) preserves ±∞, depending on det(M); Möbius maps both to the finite d/b.
function mobius(M::AbstractMatrix, above_pole::Union{Bool,Nothing}, x::T) where {T<:Rational}
    a, b, c, d = M[1,1], M[1,2], M[2,1], M[2,2]
    detM = (a*d - b*c)

    if isinf(x)
        return b == 0 ? detM * x : T(d, b)
    end

    denom = a + b * x

    if iszero(denom)
        @assert above_pole != nothing
        #If detM = 1, then above_pole -> -infty
        return above_pole ? T(-detM, 0) : T(detM, 0)
    end
    return (c + d * x) // denom
end

# Apply a per-coordinate change of basis to the region (Elower, Eupper).
# transforms[i] = [a b; c d] with det=1; acts on slope x as f_i(x) = (c + d*x)/(a + b*x).
# The pole of f_i is -a/b (absent when b=0). Splits into up to 2^n pieces per sign pattern
# of each coordinate relative to its pole. Each quadrant is first restricted via clamp, then
# the Möbius transform is applied coordinate-wise.
function basis_change(Elower::Envelope{Lower,T,D}, Eupper::Envelope{Upper,T,D},
                      transforms::Vector{M}) where {T<:Rational, D, M<:AbstractMatrix}
    n = length(transforms)

    determinants = []
    for i in 1:n
        a, b, c, d = transforms[i][1,1], transforms[i][1,2], transforms[i][2,1], transforms[i][2,2]
        @assert abs(a*d - b*c) == 1 "transform $i has determinant $(a*d - b*c), expected +-1"
        push!(determinants, a*d-b*c)
    end
    determinants = unique!(determinants)
    @assert length(unique(determinants))==1 "Can't deal with mixed determinants"
    detM = determinants[1]

    poles = [transforms[i][1,2] == 0 ? nothing : -transforms[i][1,1]//transforms[i][1,2]
             for i in 1:n]

    # Choices per coordinate: nothing means no split, true/false = above/below pole
    choices = [poles[i] === nothing ? [nothing] : [true, false] for i in 1:n]

    result = Tuple{Envelope{Lower,T,D}, Envelope{Upper,T,D}}[]

    for pattern in Iterators.product(choices...)
        # Pre-image box for this quadrant in original coordinates.
        # Above-pole (true): [pole, +∞);  below-pole (false): (-∞, pole].
        lo = T[poles[i] === nothing ? T(-1,0) : (pattern[i] ? poles[i] : T(-1,0)) for i in 1:n]
        hi = T[poles[i] === nothing ? T(1,0)  : (pattern[i] ? T(1,0)  : poles[i]) for i in 1:n]

        Elower_q, Eupper_q = clamp(Elower, Eupper, lo, hi)
        isempty(Elower_q.A) && isempty(Eupper_q.A) && continue

        transform_pt(s) = T[mobius(transforms[i], pattern[i], s[i]) for i in 1:n]

        lower_pts = [(transform_pt(s), c) for (s, c) in Elower_q.A]
        upper_pts = [(transform_pt(s), c) for (s, c) in Eupper_q.A]

        if detM==1
            push!(result, (Envelope{Lower}(lower_pts), Envelope{Upper}(upper_pts)))
        else
            push!(result, (Envelope{Lower}(upper_pts), Envelope{Upper}(lower_pts)))
        end
    end

    return result
end

# Restrict the region between two envelopes to the box [lo[i], hi[i]] for each coordinate.
# Each point is clamped to the box; push! recomputes the Pareto front.
function clamp(Elower::Envelope{Lower,T,D}, Eupper::Envelope{Upper,T,D},
               lo::Vector{T}, hi::Vector{T}) where {T, D}
    clamp_pt(p) = T[clamp(p[i], lo[i], hi[i]) for i in 1:length(p)]
    Elower_r = Envelope{Lower,T,D}()
    Eupper_r = Envelope{Upper,T,D}()
    for (p, c) in Elower.A
        push!(Elower_r, (clamp_pt(p), c))
    end
    for (p, c) in Eupper.A
        push!(Eupper_r, (clamp_pt(p), c))
    end
    return Elower_r, Eupper_r
end

function slice(v::Vector{T}, fillings::Vector{Tuple{Int,Int}}) where {T<:Real}
    return T[v[i] for i in 1:length(v) if fillings[i] == (0,0)]
end

#return a slice of an envelope with given filling slopes
function slice(E::Envelope{S,T,D}, fillings::Vector{Tuple{Int,Int}}) where {S, T<:Real, D}
    #@show fillings
    Eslice = Envelope{S,T,D}()
    for (s, c) in E.A
        if all(fillings[i] == (0,0) || strict_comp(S, fillings[i][2]//fillings[i][1], s[i]) for i in 1:length(fillings))
            push!(Eslice, (slice(s, fillings),c))
        end
    end
    return Eslice
end


#Given a lower and an upper 2D envelope, return a collection of rectangles which
#cover the region between them.
function rectangles(Elower::Envelope{Lower}, Eupper::Envelope{Upper})
    
    #sweep from left to right
    lower_pts = sort([x for (x,c) in Elower.A], by=(x -> x[1]))
    upper_pts = sort([x for (x,c) in Eupper.A], by=(x -> x[1]))
    
    ret = []
    lower_pt_queue = []

    while true
        if !isempty(upper_pts) && (isempty(lower_pts) || lower_pts[1][1] >= upper_pts[1][1])
            up = popfirst!(upper_pts)

            floor = isempty(lower_pt_queue) ? 1//0 : lower_pt_queue[end][2]


            #finish rectangles, if valid. All these points surely have lesser x coordinate
            #We need to check whether they have lesser y coordinate
            #Also we need to ensure the rectangles do not overlap.
            ceil = up[2]
            for lp in lower_pt_queue # should be sorted by x
                if all(lp .< [up[1], ceil])
                    push!(ret, (lp, [up[1], ceil]))
                    ceil = lp[2]
                end
            end

            lower_pt_queue = [[up[1], floor]]
        elseif !isempty(lower_pts)
            push!(lower_pt_queue, popfirst!(lower_pts))
        else
            break
        end
    end
    return ret
end

#Given a lower and an upper envelope, return a collection of disjoint cuboids which
#cover the region between them. Uses random pivot splits and recursion.
function cuboids(Elower::Envelope{Lower}, Eupper::Envelope{Upper})
    lower_pts = [x for (x,_) in Elower.A]
    upper_pts = [x for (x,_) in Eupper.A]

    # Prune: remove upper points that don't dominate any lower point, and vice versa
    filter!(u -> any(all(l .< u) for l in lower_pts), upper_pts)
    filter!(l -> any(all(l .< u) for u in upper_pts), lower_pts)

    isempty(lower_pts) && return Tuple{Vector,Vector}[]
    isempty(upper_pts) && return Tuple{Vector,Vector}[]

    if length(lower_pts) == 1 && length(upper_pts) == 1
        l, u = lower_pts[1], upper_pts[1]
        return all(l .< u) ? [(l, u)] : Tuple{Vector,Vector}[]
    end

    all_pts = vcat(lower_pts, upper_pts)
    pivot = rand(all_pts)
    d = length(pivot)
    k = rand(1:d)
    c = pivot[k]

    T = eltype(eltype(lower_pts))

    # Below: x_k ≤ c. Clamp upper, keep lower.
    Eupper_lo = Envelope{Upper,T,Nothing}()
    for u in upper_pts
        u_ = copy(u); u_[k] = min(u[k], c)
        push!(Eupper_lo, (u_, nothing))
    end
    Elower_lo = Envelope{Lower,T,Nothing}()
    for l in lower_pts
        push!(Elower_lo, (l, nothing))
    end

    # Above: x_k ≥ c. Clamp lower, keep upper.
    Elower_hi = Envelope{Lower,T,Nothing}()
    for l in lower_pts
        l_ = copy(l); l_[k] = max(l[k], c)
        push!(Elower_hi, (l_, nothing))
    end
    Eupper_hi = Envelope{Upper,T,Nothing}()
    for u in upper_pts
        push!(Eupper_hi, (u, nothing))
    end

    return vcat(cuboids(Elower_lo, Eupper_lo), cuboids(Elower_hi, Eupper_hi))
end


include("BasisChange.jl")

end
