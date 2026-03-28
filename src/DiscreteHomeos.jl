@enum RoundMode DOWN=1 UP=2
@enum DIR LEFT=1 RIGHT=2

import Base: *

getindex(t::Tuple, d::DIR) = (d==LEFT ? t[1] : t[2])

struct DiscreteHomeo{T} <: Homeo #todo: precompute the output heights
    ordering::Vector{Tuple{T,T}}
    dir::DIR
    roundmode::RoundMode
end

struct DiscreteHomeoAlt{T,S} <: Homeo
    ordering_l::Vector{Tuple{T,S}}
    ordering_r::Vector{Tuple{T,S}} #T is the type of the inputs, S is the type of the middle index
    roundmode::RoundMode
end

function _is_valid(f::DiscreteHomeoAlt)
    ol = [x[1] for x in f.ordering_l]
    ol_m = [x[2] for x in f.ordering_l]
    or = [x[1] for x in f.ordering_r]
    or_m = [x[2] for x in f.ordering_r]
    return issorted(ol) && issorted(ol_m) && issorted(or) && issorted(or_m) &&
    length(unique(ol))==length(ol) &&
    length(unique(or))==length(or) &&
    unique(ol_m) == unique(or_m)
end

function fix_heights!(c::Cand{H}) where {S,T,H<: DiscreteHomeoAlt{T,S}}
    for fnum in 0:2*c.bt.ntets-1
        #gather all the heights
        heights = Rational{Int}[]
        
        for k in 0:2
            i, J = c.bt.forward[(fnum, k)]

            f = c[J]
            for tup in filter(tup -> tup[1][1]==i, f.ordering_l)
                push!(heights, tup[1][2])
            end
        end

        heights = sort(unique(heights))
        for k in 0:2
            i, J = c.bt.forward[(fnum, k)]
            for h in heights
                c[J] = insert_element(c[J], (i,h))
            end

            for h in heights
                @assert length(
                searchsorted(c[J].ordering_l, ((i,h),(i,h)), by=x->x[1])
               )>=1
            end
        end
    end
end

function add_height!(c::Cand{H}, fnum::Int, height) where {H}
    for k in 0:2
        i, J = c.bt.forward[(fnum, k)]
        c[J] = insert_element(c[J], (i,height))
    end
end


function insert_element(f::DiscreteHomeoAlt{T,S},el) where {T,S}
    i = searchsortedlast(f.ordering_l, (el, zero(S)), by=x->x[1])
    @assert i >= 1
    if f.ordering_l[i][1] == el
        return f
    else
        ordering_l = copy(f.ordering_l)
        #find the closest element 
        push!(ordering_l, (el, f.ordering_l[i][2]))
        sort!(ordering_l)
        return DiscreteHomeoAlt(ordering_l, f.ordering_r, f.roundmode)
    end
end

function inv(f::DiscreteHomeoAlt)
    return DiscreteHomeoAlt(f.ordering_r, f.ordering_l, f.roundmode)
end

struct CompositeHomeo{H} <: Homeo
    l::H
    r::H
end

function *(h1::Homeo, h2::Homeo)
    return CompositeHomeo(h1,h2)
end

function (f::CompositeHomeo)(y)
    return f.r(f.l(y))
end

function inv(c::CompositeHomeo)
    return CompositeHomeo(inv(c.r), inv(c.l))
end

function inv(d::DIR)
    if d==LEFT
        return RIGHT
    else
        return LEFT
    end
end

function apply_dir(d::DIR, tup::Tuple{T,T}) where {T}
    if d == LEFT
        return tup
    else
        return (tup[2],tup[1])
    end
end

function (d::DIR)(tup::Tuple{T,T}) where {T}
    return apply_dir(d,tup)
end


myzero(T::Type{Tuple{A,B}}) where {A,B} = (myzero(A),one(B))
myzero(X::Type{T}) where T<:Real= zero(X)
myzero(::Type{Rational{T}}) where {T} = zero(Rational{T})

function is_low(f::DiscreteHomeoAlt{T,S}, x::T) where {T,S}
    r=searchsorted(f.ordering_l, (x,myzero(S)), by = a->a[1])
    @assert length(r)==1


    mid_index = f.ordering_l[r[1]][2]

    s=searchsorted(f.ordering_l,  (myzero(T), mid_index), by = a->a[2])
    return r[1]==s[1]
end

function (f::DiscreteHomeoAlt{T,S})(x::T) where {T,S}
    r=searchsorted(f.ordering_l, (x,myzero(S)), by = a->a[1])
    @assert length(r)==1

    mid_index = f.ordering_l[r[1]][2]
    s=searchsorted(f.ordering_r,  (myzero(T), mid_index), by = a->a[2])
    @assert length(s) >= 1

    if f.roundmode==DOWN
        return f.ordering_r[s[1]][1]
    else
        return f.ordering_r[s[end]][1]
    end
end

function (f::DiscreteHomeo{T})(x::T) where {T}
    if length(f.ordering)==0
        return myzero(T)
    end

    if f.roundmode == DOWN
        r = searchsortedfirst(f.ordering, (x,x), by=a->a[f.dir])
        if 1<=r<=length(f.ordering) && f.ordering[r][f.dir] != x
            r=r-1
        end
        return f.ordering[min(max(r,1),length(f.ordering))][inv(f.dir)]
    else
        r = searchsortedlast(f.ordering, (x,x), by=a->a[f.dir])
        if 1<=r<=length(f.ordering) && f.ordering[r][f.dir] != x
            r=r+1
        end
        return f.ordering[min(max(r,1),length(f.ordering))][inv(f.dir)]
    end
end


#=
function insert_left!(f::DiscreteHomeo{T},el) where {T}
    r=searchsortedfirst(f.ordering, (el,el), by=a->a[f.dir])
    t = f.ordering[r]
    s = f.ordering[r+1]
    tp = (el, t[inv(f.dir)])
    sp = (el, s[inv(f.dir)])
    insert!(f.ordering(f.ordering, r, tp))
    insert!(f.ordering(f.ordering, r+1, sp))
end

function insert_right!(f::DiscreteHomeo{T},el) where {T}
    return inv(insert_left!(inv(f),el))
end



function splittings(f::DiscreteHomeo{T}, i::T, i1::T, i2::T) where {T}
    o=unique(f.ordering)
    index_range = searchsorted(o, (i,i), by=x->x[f.dir])

    ret = []

    @assert length(index_range) >= 0

    for ind in index_range
        t=copy(o) #makes a copy

        for j in index_range
            t[j] = apply_dir(f.dir, (j <= ind ? i1 : i2, t[j][inv(dir)]) )
        end
        insert!(t, ind+1, apply_dir(f.dir, (i2, t[ind][inv(f.dir)])))

        push!(ret, cleanup(DiscreteHomeo(t, f.dir, f.roundmode)))
    end
    return ret
end
=#

#=
function bump_left(ord::Vector, i, h)
    o = copy(ord)
    for i in eachindex(o)
        o[i] = f.dir( (bump(i, h, o[i][f.dir]), o[i][inv(f.dir)]) )
    end 
    return o
end
function bump_right(ord::Vector, i, h)
    o=copy(ord)
    for i in eachindex(o)
        o[i] = f.dir( (o[i][f.dir], bump(i, h, o[i][inv(f.dir)])) )
    end 
    return o
end

function bump_left(f::DiscreteHomeo, i, h)
    o=copy(d.ordering)
    for i in eachindex(o)
        o[i] = f.dir( (bump(i, h, o[i][f.dir]), o[i][inv(f.dir)]) )
    end 

    return DiscreteHomeo(o, f.dir, f.roundmode)
end
function bump_right(f::DiscreteHomeo, i, h)
    o=copy(d.ordering)
    for i in eachindex(o)
        o[i] = f.dir( (o[i][f.dir], bump(i, h, o[i][inv(f.dir)])) )
    end 

    return DiscreteHomeo(o, f.dir, f.roundmode)
end

function bump(i::Int, h::Int, el::Tuple{Int,Int})
    if el[1] == i
        if el[2] <= h
            return el
        else
            return (i, el[2] + 1)
        end
    end
end
=#

function midpoint(i::Rational{Int}, j::Union{Rational{Int},Nothing})
    if j==nothing
        return i+1
    else
        return (i+j)//2
    end
end

function midpoint(i::Tuple{Int, Rational{Int}}, j::Union{Tuple{Int,Rational{Int}},Nothing})
    if j==nothing || i[1] != j[1]
        (i[1], i[2]+1)
    else
        @assert i <= j
        (i[1], (i[2] + j[2])//2)
    end
end

function midheights(f::DiscreteHomeoAlt{T,S}) where {T,S}
    return unique([x[2] for x in f.ordering_l])
end

function splittings(c::Cand{H}, s::State) where {S,T,H<: DiscreteHomeoAlt{T,S}}

    h = s.x #height inside the track
    t = s.e

    _i, J = c.bt.forward[t]
    f=c[J]

    search_l = searchsorted(f.ordering_l, ((_i,h),myzero(S)), by=x->x[1])
    @assert length(search_l) ==1
    ind = search_l[1]
    
    _, height_mid = f.ordering_l[ind] 
    index_range_l = searchsorted(f.ordering_l, (myzero(T),height_mid), by=x->x[2])
    index_range_r = searchsorted(f.ordering_r, (myzero(T),height_mid), by=x->x[2])

    if ind==index_range_l[1] #no splitting required
        @assert is_low(f, (_i,h))
        return [c]
    else
        @assert !is_low(f, (_i,h))
    end

    height_mid_next = (index_range_r.stop < length(f.ordering_r) ? f.ordering_r[index_range_r.stop+1][2] : nothing)
    height_mid_new = midpoint(height_mid, height_mid_next)
    @assert height_mid < height_mid_new
    @assert height_mid_new != height_mid_next

    ret = []
    for i in index_range_r
        cnew = copy(c)
        next_index = midpoint(f.ordering_r[i][1], i < length(f.ordering_r) ? f.ordering_r[i+1][1] : nothing)

        ordering_l = copy(f.ordering_l)
        ordering_r = copy(f.ordering_r)


        for j in ind:index_range_l.stop
            @assert ordering_l[j][2] == height_mid
           ordering_l[j] = (ordering_l[j][1], height_mid_new)
        end
        @assert issorted(ordering_l)

        @assert unique([x[2] for x in ordering_l]) == sort(unique(vcat([x[2] for x in f.ordering_l], [height_mid_new])))

        for j in i+1:index_range_r.stop
            @assert ordering_r[j][2] == height_mid
            ordering_r[j] = (ordering_r[j][1], height_mid_new)
        end
        @assert issorted(ordering_r)

        push!(ordering_r, (next_index, height_mid_new))
        sort!(ordering_r)

        @assert unique([x[2] for x in ordering_r]) == sort(unique(vcat([x[2] for x in f.ordering_r], [height_mid_new])))

        fnew = DiscreteHomeoAlt(ordering_l, ordering_r, f.roundmode)

        @assert _is_valid(fnew)

        cnew[J] = fnew
        #fix_heights!(cnew)

        #find the n
        add_height!(cnew, cnew.bt.backwardfan[(next_index[1], J)][1], next_index[2])

        @assert is_low(cnew[J], (_i,h))
        push!(ret, cnew)
    end

    return ret
end

#=
function splittings(c::Cand{H}, s::State) where {H <: CompositeHomeo}
    #do this in two steps. First we add a new internal height to the nearby junction. Then we split apart on the other side in one of several ways
    #
    h = s.x #height inside the track
    t = s.e

    i, J = c.bt.forward[t]
    f=c[J]
    l=f.l
    r=f.r

    #now let's try to extend


    o = unique(l.ordering)
    index_range_l = searchsorted(o, ((i,h),(i,h)), by=x->x[l.dir])
    @assert length(index_range_l) == 1 #this is because f.l should be compressing
    ind = index_range_l[1]
    y = o[ind][inv(l.dir)]

    index_range_mid = searchsorted(o, (y,y), by=x->x[inv(l.dir)])
    @assert ind in index_range_mid
    if ind == index_range_mid[1] #if this is the lowest in the range, then it already passes through
        return [c]
    end
    #now we know which range we want to split in index_range_mid
    #


    lnew = cleanup(DiscreteHomeo(o_bumped, l.dir, l.roundmode))
    rnew = bump_left(r, y[1], y[2])

    index_rand_mid_1 = o[l.ordering]

    o_r = unique(r.ordering)
    index_range_r = searchsorted(o_r, (y,y), by=x->x[r.dir])
    @assert length(index_range_r) >= 1 #these are all the possible choices for the splitting

    for indr in index_range_r
        c2 = copy(c)
        l = bump_right()

        bump()
        #need to shift up y in both l and r
        #need to split apart indr
    end
end

function is_contracting(f::DiscreteHomeo)
    return length(unique(x[f.dir] for x in f.ordering)) == length(unique(f.ordering))
end
function is_expanding(f::DiscreteHomeo)
    return length(unique(x[inv(f.dir)] for x in f.ordering)) == length(unique(f.ordering))
end

function split!(c::Cand{H}, fnum, h)  where {H<:CompositeHomeo}
    @assert height in c.heights[fnum]
    #do the relabel.
    relabel = Dict(i=> (i <= h ? i : i+1) for i in c.heights[fnum])
    c.heights[fnum] = [relabel[i] for i in c.heights[fnum]]
    push!(c.heights[fnum], h+1)
    sort!(c.heights[fnum])

    for k in 0:2
        i, J = c.bt.forward[(fnum, k)]
        f = copy(c[J].l)

        for j in eachindex(f.ordering)
            x,y = f.dir(f.ordering[j])
            #relabel (i, y) to (i,y+1) whenever y >= h
            f.ordering[j] = f.dir(((x[1], (x[1] != i || x[2] <= h) ? x[2] : x[2] + 1) , y))
        end

        #now we have to insert (i, h+1) somewhere.
        index_range = searchsorted(f.o, ((i,h),(i,h)), by=x->x[f.dir])
        @assert length(index_range)==1
        ind = index_range[1]

        x,y = f.dir(f.o[ind])
        @assert x == (i,h)
        push!(f.o, f.dir((i,h+1),y))
        c[J] = cleanup(f) * c[J].r  
    end
    return c
end

function splittings(c::Cand, fnum::Int, h::Int)
    c2 = copy(c)
    relabel_heights!(c, fnum, Dict(i=> (i<= h ? i : i+1) for i in c.heights[fnum]) )
    return splittings(c2, fnum, h, h, h+1)
end

function splittings(c::Cand{DiscreteHomeo{T}}, fnum::Int, h::Int, h1::Int, h2::Int) where {T}
    ret = Cand{DiscreteHomeo{T}}[]
    i0, J0 = c.bt.forward[(fnum,0)]
    i1, J1 = c.bt.forward[(fnum,1)]
    i2, J2 = c.bt.forward[(fnum,2)]

    #we need to find all the possible splittings in these three junctions

    for f0 in splittings(c[J0], (i0,h), (i0,h1), (i0,h2))
        c0 = copy(c)
        c0[J0] = f0
        for f1 in splittings(c0[J1], (i1, h), (i1,h1), (i1, h2))
            c1 = copy(c0)
            c1[J1] = f1
            for f2 in splittings(c2[J2], (i2,h), (i2,h1), (i2,h2))
                c2 = copy(c1)
                c2[J2]=f2
                push!(ret, c2)
            end
        end
    end
    return ret
end

function relabel_heights!(c::Cand{DiscreteHomeo{T}},fnum::Int, relabel::Dict{Int,Int}) where {T}
    c.heights[fnum] = [relabel[i] for i in c.heights[fnum]]

    for i in 0:2
        tr = (fnum, i)
        j, J = c.bt.forward[tr]
        f = c[J]
        for i in eachindex(f.ordering)
            x = f.ordering[i][f.dir]
            y = f.ordering[i][inv(f.dir)]
            if x[1] == j #if this edge is coming in from tr
                f.ordering[i] = apply_dir(f.dir, ((j,relabel[x[2]]),y))
            end
        end
    end
end
=#


function inv(f::DiscreteHomeo)
    return DiscreteHomeo(f.ordering, inv(f.dir), f.roundmode)
end


function jiggle(f::DiscreteHomeo, r::T) where {T}
    ordering = copy(f.ordering)
    #@assert issorted(ordering)

    for i in 2:length(ordering)-1
        if rand() < r &&
           xor(ordering[i-1][1] == ordering[i][1],ordering[i][1] == ordering[i+1][1]) &&
           xor(ordering[i-1][2] == ordering[i][2],ordering[i][2] == ordering[i+1][2]) &&
           ordering[i] != ordering[i+1]

            s=rand()
            if s < 1/3
                ordering[i] = ordering[i-1]
            elseif s < 2/3
                ordering[i] = (ordering[i-1][1], ordering[i+1][2])
            else
                ordering[i] = (ordering[i+1][1], ordering[i-1][2])
            end
        end
    end
    #@assert issorted(ordering)
    #@assert unique([x[1] for x in ordering]) == unique([x[1] for x in f.ordering])
    #@assert unique([x[2] for x in ordering]) == unique([x[2] for x in f.ordering])

    #delete duplicates
    
    return DiscreteHomeo(ordering, f.dir, f.roundmode)
end


function longitude_to_candidate(bt::BoundaryTriangulation, l)
    d=Dict{Junction, DiscreteHomeo{Tuple{Int,Int}}}()
    for J in bt.junctions 
        left = Tuple{Int,Int}[]
        right = Tuple{Int,Int}[]

        for k in 0:J.left_len-1
            tri_num, index = bt.forwardfan[(k,J)] :: Track
            push!(left,(k,0))
            for j in 1:l[tri_num]
                push!(left, (k,j))
            end
        end

        for k in 0:J.right_len-1
            tri_num, index = bt.backwardfan[(k,J)] :: Track
            push!(right,(k,0))
            for j in 1:l[tri_num]
                push!(right,(k,j))
            end
        end

        ordering = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
        push!(ordering, (popfirst!(left),popfirst!(right)))

        #annoying part: we need to add some dummy edges that fill up the edges which have coefficient 0

        while length(left) > 0 || length(right) > 0
            if length(left) > 0 && left[1][2] == 0
                push!(ordering, (popfirst!(left), ordering[end][2]))
            elseif length(right) > 0 && right[1][2] == 0
                push!(ordering, (ordering[end][1],popfirst!(right)))
            else
                @assert length(left) > 0 && length(right) > 0
                push!(ordering, (popfirst!(left), popfirst!(right)))
            end
        end
        

        H = cleanup(DiscreteHomeo(ordering, LEFT, DOWN))
        d[J] = H
        d[inv(J)] = inv(d[J])
    end

    return Cand(bt, ArrayDict(d)) #make into offset array
end


function random_ordering(left::Vector{T}, right::Vector{T}) where {T}
    @assert length(left)>0
    @assert length(right)>0
    if length(left)==1
        return Tuple{T,T}[(left[1],r) for r in right]
    elseif length(right)==1
        return Tuple{T,T}[(l,right[1]) for l in left]
    else
        k1 = rand(1:length(left)-1)
        k2 = rand(1:length(right)-1)
        return vcat(random_ordering(left[1:k1], right[1:k2]), random_ordering(left[k1+1:end], right[k2+1:end]))
    end
end

function all_orderings(left::Vector{T}, right::Vector{T}) where {T}
    @assert issorted(left)
    @assert issorted(right)
    @assert length(left) > 0
    @assert length(right) > 0
    if length(left)==1
        return [Tuple{T,T}[(left[1],r) for r in right]]
    elseif length(right)==1
        return [Tuple{T,T}[(l,right[1]) for l in left]]
    else
        ret = []
        for x in all_orderings(left[2:end], right)
            push!(ret, vcat(Tuple{T,T}[(left[1], right[1])], x))
        end
        for x in all_orderings(left, right[2:end])
            push!(ret, vcat(Tuple{T,T}[(left[1], right[1])], x))
        end
        return ret
    end
end

function lrheights(bt::BoundaryTriangulation, J::Junction, heights)
    lheights = Tuple{Int,Int}[]
    rheights = Tuple{Int,Int}[]

    for i in 0:J.left_len-1
        track = bt.forwardfan[(i,J)]
        index, _ = track
        for h in heights[index]
            push!(lheights, (i, h))
        end
    end
    for i in 0:J.right_len-1
        track = bt.backwardfan[(i,J)]
        index, _ = track
        for h in heights[index]
            push!(rheights, (i, h))
        end
    end
    return lheights, rheights
end

function all_cands(bt, roundmode)
    heights = OffsetArrays.Origin(0)([[0] for i in 1:2*bt.ntets])

    homeolists = [
                  [cleanup(DiscreteHomeo{Tuple{Int,Int}}(o,LEFT,roundmode))
                   for o in all_orderings(lrheights(bt, J, heights)...)]
                   for J in bt.junctions]


    println("ncands: $(prod(map(length, homeolists)))")

    gen = (
           begin
        d=Dict{Junction, DiscreteHomeo{Tuple{Int,Int}}}()
            for (J,f) in zip(bt.junctions,l)
                d[J] = f
                d[inv(J)] = inv(f)
            end 
            Cand(bt, ArrayDict(d))
           end
    for l in Iterators.product(homeolists...)
   )

    return gen
end

function random_discrete_homeo(J::Junction, thickness::Int, roundmode::RoundMode)
    ret = cleanup(DiscreteHomeo{Tuple{Int,Int}}(
                                 random_ordering(
                                                 sort!(collect((i,j) for i in 0:J.left_len-1, j in 1:thickness)[1:end]),
                                                 sort!(collect((i,j) for i in 0:J.right_len-1, j in 1:thickness)[1:end])),
                                 LEFT, roundmode)) 
end

function basic_homeo(J::Junction, roundmode::RoundMode)
    left = Tuple{Tuple{Int,Rational{Int}},Int}[((i,0), 0) for i in 0:J.left_len-1]
    right = Tuple{Tuple{Int,Rational{Int}},Int}[((i,0), 0) for i in 0:J.right_len-1]

    return DiscreteHomeoAlt{Tuple{Int,Rational{Int}},Rational{Int}}(left, right, roundmode)
end


function trace_forwards(s::State{T}, c::Cand) where {T}
	i1,J = c.bt.forward[s.e]
	
    (i2,x2) = c[J]((i1,s.x))

    return State{T}(x2, c.bt.backwardfan[(i2,J)])
end

function verify_low(s::State{T}, c::Cand, n::Int) where {T}
    #@show "verifying $(n)"
    for i in 1:n
        i1, J = c.bt.forward[s.e]
        @assert is_low(c[J], (i1,s.x))
        (i2,x2) = c[J]((i1,s.x))
        s = State{T}(x2, c.bt.backwardfan[(i2,J)])
    end
end

function assign_height_left(f::DiscreteHomeoAlt{T}, r::T; roundmode=:mean) where {T}
    ind_range = searchsorted(f.ordering_l, (r,r), by = x->x[1])
    @assert length(ind_range) == 1
    ind = ind_range[1]
    mid_height = f.ordering_l[ind][2]

    #how many heights are below this guy?



    l2 = unique([x[2] for x in f.ordering_l])
    l1 = collect(filter(x->x<= mid_height, l2))

    x2 = collect(filter(x->x[2] == mid_height, f.ordering_l))
    x1 = collect(filter(x->x[1] <= r, x2))

    low = (length(l1)-1)/length(l2) + (length(x1)-1)/(length(x2)*length(l2))
    high = (length(l1)-1)/length(l2) + (length(x1))/(length(x2)*length(l2))

    if roundmode==:mean
        return (low+high)/2
    elseif roundmode==:up
        return high
    elseif roundmode==:down
        return low
    end
    @assert false
end

function assign_height_right(f::DiscreteHomeoAlt{T}, r::T; roundmode=:mean) where {T}
    return assign_height_left(inv(f), r; roundmode=roundmode)
end

function assign_height_left(f::DiscreteHomeo{T}, r::T; roundmode=:mean) where {T}
    o=unique(f.ordering)
    index_range = searchsorted(o, (r,r), by=x->x[f.dir])
    if roundmode==:mean
        (index_range.start + index_range.stop-1)/(2*length(o))
    elseif roundmode==:up
        (index_range.stop)/length(o)
    elseif roundmode==:down
        (index_range.start-1)/length(o)
    end
end

function assign_height_right(f::DiscreteHomeo{T}, r::T; roundmode=:mean) where {T}
    o=unique(f.ordering)
    index_range = searchsorted(o, (r,r), by=x->x[inv(f.dir)])
    if roundmode==:mean
        (index_range.start + index_range.stop-1)/(2*length(o))
    elseif roundmode==:up
        (index_range.stop)/length(o)
    elseif roundmode==:down
        (index_range.start-1)/length(o)
    end
end

function cleanup(f::DiscreteHomeo{T}) where {T}
    o = Tuple{T,T}[]
    for i in unique(f.ordering)
        if length(o) > 1 && o[end][1] != i[1] && o[end][2] != i[2]
            push!(o,o[end])
        end
        push!(o,i)
    end
    return DiscreteHomeo(o, f.dir, f.roundmode)
end


function basic_cand(bt, roundmode)
    d=Dict{Junction, DiscreteHomeoAlt{Tuple{Int,Rational{Int}}, Rational{Int}}}()

    for j in bt.junctions
        d[j] = basic_homeo(j, roundmode)
        d[inv(j)] = inv(d[j])
    end

    return Cand(bt, ArrayDict(d))

end

function heights(c::Cand{X}, fnum::Int) where {X<:DiscreteHomeoAlt}
    ret = Rational{Int}[]
    for k in 0:2
        i, J = c.bt.forward[(fnum, k)]

        f = c[J]
        for tup in filter(tup -> tup[1][1]==i, f.ordering_l)
            push!(ret, tup[1][2])
        end
    end
    return unique(ret)
end
function heights(c::Cand{X}, fnum::Int) where {X<:DiscreteHomeo}
    ret = Int[]
    for k in 0:2
        i, J = c.bt.forward[(fnum, k)]

        f = c[J]
        for tup in filter(tup -> tup[f.dir][1]==i, f.ordering)
            push!(ret, tup[f.dir][2])
        end
    end
    return unique(ret)
end
