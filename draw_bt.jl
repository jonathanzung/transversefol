import Luxor
using Subscripts

MyPoint = SVector{2,Rational{Int}}
MyFloatPoint = SVector{2,Float64}


function coboundary(bt, A::DefaultDict{Junction, T}) where {T}
    ret = DefaultDict{Track,T}(zero(T))

    for (J,weight) in A
        for t in right_tracks(bt,J)
            ret[t] += weight
        end
        for t in left_tracks(bt,J)
            ret[t] -= weight
        end
    end
    return ret
end

function bdry(bt, A::DefaultDict{Track, T}; scale=x->one(T)) where {T}
    ret = DefaultDict{Junction,T}(zero(T))
    for (tr,weight) in A
        _,J = bt.forward[tr]
        ret[J] += scale(tr)*weight
        _,J = bt.backward[tr]
        ret[J] -= scale(tr)*weight
    end
    return ret
end

function draw(c::Cand{Hm}, i::Int; SCALE=200, PAD=200, JUNCTIONHEIGHT=1//8, curve=nothing) where {Hm} #i is which cusp to draw
    #curve is an ArrayDict giving weights of a curve to draw on the boundary.
    #weight should be a map from Tracks to displacement vectors
    #
    
    bt = c.bt
    poles = ladderpoles(c.bt,i)

    H = maximum(map(length, poles))
    W = length(poles)
    Point=Luxor.Point

    SCALEPT = MyPoint(W*SCALE, H*SCALE)

    offset = DefaultDict{Junction, MyPoint}(MyPoint(0,0))

    pole_index = Dict{Junction, Int}()
    pole_edges = Set{Track}()

    for (j, pole) in enumerate(poles)
        #pole is a sequence of edges

        for (i,track) in enumerate(pole)
            _,J = bt.backward[track]
            offset[J] = MyPoint((mod(-j+1,length(poles)))//length(poles), 0)
            pole_index[J]=j
            push!(pole_edges, track)
        end
    end

    displacement = coboundary(bt, offset)

    for track in keys(displacement)
        displacement[track] += bt.weights[track]
    end

    for i in 1:50
        #we can't round displacements, or else it might not be an coexact cocycle anymore
        #but it's okay to round the forces.

        springforce = bdry(bt,displacement; scale = x-> x in pole_edges ? 1 : 1//3)
        for J in keys(springforce)
            springforce[J] = round.(Int, 1024 .* springforce[J])//1024
        end

        delta = coboundary(bt,springforce)

        for track in keys(displacement)
            displacement[track] += [0,1//8] .* delta[track]
        end
    end

    Luxor.Drawing(SCALEPT[1]+2*PAD, SCALEPT[2]+2*PAD, "cusp$(i).svg")
    Luxor.origin(PAD,PAD)
    #Luxor.background("white")
    Luxor.box(-Point(PAD-1, PAD-1), Point((SCALEPT .+ (PAD - 1))...), action = :clip)
    Luxor.sethue("black")
    Luxor.setline(0.25)


    function control(J::Junction)
        if isodd(pole_index[J])
            SCALE*MyPoint(1//4,1//4)
        else
            SCALE*MyPoint(1//4,-1//4)
        end
    end

    function vert_control(J::Junction)
        if isodd(pole_index[J])
            SCALE*MyPoint(-JUNCTIONHEIGHT,JUNCTIONHEIGHT)
        else
            SCALE*MyPoint(JUNCTIONHEIGHT,JUNCTIONHEIGHT)
        end

    end

    function draw_edge(tr,basepoint) #draw a track, starting from a given basepoint
        #=
        i1,J1 = bt.backward[tr]
        i2,J2 = bt.forward[tr]
        left_heights = collect(filter(y->y[1] == i1, collect(x[inv(c[J1].dir)] for x in c[J1].ordering))) 
        right_heights = collect(filter(y->y[1] == i2, collect(x[c[J2].dir] for x in c[J2].ordering))) 
        for h in unique([h[2] for h in vcat(left_heights, right_heights)])
        =#

        if curve!=nothing
            Luxor.setline(1+2*abs(curve[tr][1]))
            draw_edge(tr,  basepoint, 0; roundmode=:mean)
            Luxor.setline(2)
        else
            for h in heights(c, tr[1])
                draw_edge(tr, basepoint, h;roundmode=:up)
                draw_edge(tr, basepoint, h;roundmode=:down)
            end
        end
    end

    function draw_edge(tr,basepoint, height; roundmode=:mean)
        i1,J1 = bt.backward[tr]
        i2,J2 = bt.forward[tr]

        p1 = SCALEPT.*(basepoint) + (assign_height_right(c[J1], (i1,height); roundmode=roundmode) - 0.5)*vert_control(J1)
        p2 = p1 .+ control(J1)

        p4 = SCALEPT.*(basepoint + displacement[tr]) + (assign_height_left(c[J2], (i2,height);roundmode=roundmode) - 0.5)*vert_control(J2)
        p3 = p4 .- control(J2)

        Luxor.move(p1)
        Luxor.curve(p2,p3,p4)
        Luxor.strokepath()
    end

    function draw_junction(J,basepoint)
        _pt = Point((SCALEPT.*basepoint)...)
        #circ=Luxor.circle(_pt, 5)
        Luxor.fontsize(18)
        Luxor.label("J"*sub("$(J.index)")*super(J.inv ? "-1" : ""), 
                    isodd(pole_index[J]) ? :NE : :NW, _pt, leader=false, offset=SCALE*0.1)
        #Luxor.label(L"x^2", isodd(pole_index[J]) ? :NE : :NW, _pt, leader=false, offset=SCALE*0.1)
    end

    #now let's do a dfs and draw all the edges

    _, Jstart = bt.backward[poles[1][1]]
    visited=Set{Tuple{Junction,MyPoint}}() #we need to do this with rationals
    Q=Queue{Tuple{Junction,MyPoint}}()
    enqueue!(Q, (Jstart, MyPoint(-3//4,-3//4)))
	while !isempty(Q)
		J,pt = dequeue!(Q)
        if (J,pt) in visited
            continue
        end
        push!(visited,(J,pt))
        draw_junction(J,pt)

        for track in right_tracks(bt,J)
            _,J2 = bt.forward[track]

            next = (J2, pt + displacement[track])
            draw_edge(track, pt)

            if all(abs.(next[2]) .< 2)
                enqueue!(Q, next)
            end
        end
	end

    #now let's draw some random paths

    Luxor.sethue("orange")
    Luxor.setline(2)
 
    for l in -10:10
        s=State(rand_init(Hm), poles[1][1])
        pt = MyPoint(-3//4,-3//4+l)

        for i in 1:80
            draw_edge(s.e, pt, s.x)
            pt = pt + displacement[s.e]
            s=trace_forwards(s, c)
        end
    end


    Luxor.sethue("green")
    Luxor.box(Point(0,0), Point(SCALEPT...), action = :stroke)
    Luxor.clipreset()

    Luxor.finish()
    Luxor.preview()


end
