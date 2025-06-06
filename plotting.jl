#using Plots
using PlotlyJS
#plotlyjs()
#
#

function plotjs(A::Vector{Envelope})
	PlotlyJS.plot([_plotjs(E) for E in A])
end

function plotjs(E::Envelope)
	plotjs(Envelope[E])
end
const TAUT_COLOUR = "#00cc96"
const POS_CONTACT_COLOUR = "#eeee00"
const NEG_CONTACT_COLOUR = "#00eeee"
const OBSTRUCTION_COLOUR ="#ef553b" 
const LONGITUDE_COLOUR = "#636EFA"

function _plotjs(E::Envelope{S}; color=TAUT_COLOUR, name="") where {S<:Union{Upper,Lower}}
	pts = [x[1] for x in E.A]
	dim = length(pts[1])

	if dim==2
        marker = if color==nothing
            attr()
        else
            attr(color=color)
        end
		all_pts = staircase(E)
		PlotlyJS.scatter(x=[x[1] for x in all_pts],y=[x[2] for x in all_pts], mode="lines", line=attr(color=color), name=name, legendgroup=name)
	elseif dim==3
        marker = if color==nothing
            attr(size=3)
        else
            attr(color=color,size=3)
        end
		PlotlyJS.scatter(x=[x[1] for x in pts],y=[x[2] for x in pts], z=[x[3] for x in pts], mode="markers", type="scatter3d", marker=marker)
	else
		@assert false
	end
end

function _plotjs(E::Envelope{S,T}; color=nothing) where {S,T}
	pts = filter(inbounds, [x[1] for x in E.A])


	if length(pts)==0
		return PlotlyJS.scatter(x=T[],y=T[], mode="markers")
	end
	dim = length(pts[1])

	if dim==2
        _pts = filter(x->abs(x[1]) <= CLIP && abs(x[2]) <= CLIP, pts)
        marker = if color==nothing
            attr()
        else
            attr(color=color)
        end
		PlotlyJS.scatter(x=T[x[1] for x in _pts],y=T[x[2] for x in _pts], mode="markers", marker=marker)
	elseif dim==3
        marker = if color==nothing
            attr(size=3)
        else
            attr(color=color,size=3)
        end

		PlotlyJS.scatter(x=T[x[1] for x in pts],y=T[x[2] for x in pts], z=T[x[3] for x in pts], mode="markers", type="scatter3d", marker=marker)
	else
		@assert false	
	end
end

function _plotjs(E1::Envelope{Lower}, E2::Envelope{Upper}; color=TAUT_COLOUR, name="")
	pts1 = [x[1] for x in E1.A]
	dim = length(pts1[1])

	#pts2 = [x[1] for x in E2.A]
	#=

	@assert length(pts1[1])==2

	sort!(pts1, by=x->x[1])
	sort!(pts2, by=x->x[1])
	all_pts1=sort(vcat(pts1,[(pts1[i+1][1], pts1[i][2]) for i in 1:length(pts1)-1]), by=x->(x[1],-x[2]))
	all_pts2=sort(vcat(pts2,[(pts2[i][1], pts2[i+1][2]) for i in 1:length(pts2)-1]), by=x->(x[1],-x[2]))
	=#

	if dim==2
		all_pts1=staircase(E1)
		all_pts2=staircase(E2)

        if length(all_pts1)==0
            push!(all_pts1, [CLIP,CLIP])
        end
        if length(all_pts2)==0
            push!(all_pts2, [-CLIP,-CLIP])
        end

		if all(all_pts1[1] .<= all_pts2[1])
			pushfirst!(all_pts2, [all_pts1[1][1], all_pts2[1][2]])
			pushfirst!(all_pts1, [all_pts1[1][1], all_pts2[1][2]])
		end

		if all(all_pts1[end] .<= all_pts2[end])
			push!(all_pts1, [all_pts2[end][1], all_pts1[end][2]])
			push!(all_pts2, [all_pts2[end][1], all_pts1[end][2]])
		end

		return [PlotlyJS.scatter(x=[x[1] for x in all_pts1],y=[x[2] for x in all_pts1], mode="lines", fill="tonexty", name=name, legendgroup = name, fillcolor="#00000000", line=attr(color=color)),
				PlotlyJS.scatter(x=[x[1] for x in all_pts2],y=[x[2] for x in all_pts2], mode="lines", fill="tonexty", name=name, legendgroup = name, line=attr(color=color)),
				]
	elseif dim==3
        cubes = []
        for (p1,_) in E1.A
            for (p2,_) in E2.A
                if all(p1 .<= p2)
                    push!(cubes, cube(p1,p2))
                end
            end
        end
        return [combine_meshes(cubes)]
		#return [_plotjs(E1; color=color), _plotjs(E2, color=color)]
	end
end

function combine_meshes(v::Vector) #a vector of tuples
    xs = collect(Iterators.flatten(m.x for m in v))
    ys = collect(Iterators.flatten(m.y for m in v))
    zs = collect(Iterators.flatten(m.z for m in v))
    facecolors = collect(Iterators.flatten(m.facecolor for m in v))

    is = Int[]
    js = Int[]
    ks = Int[]
    curroffset = 0

    for m in v
        append!(is, m.i .+ curroffset)
        append!(js, m.j .+ curroffset)
        append!(ks, m.k .+ curroffset)
        curroffset += length(m.x)
    end

    return mesh3d(x=xs,y=ys,z=zs,i=is,j=js,k=ks,facecolor=facecolors, flatshading=true, showlegend=true)
end

function cube(p1, p2)
    facecolor = repeat([
    	#"rgb(50, 200, 200)", #GB
    	#"rgb(100, 200, 255)", #BBG
    	"rgb(200, 200, 50)", #Yellow
    	"rgb(150, 200, 115)", #G
    	"rgb(150, 200, 115)", #G
    	#"rgb(100, 200, 255)", #BBG
    	"rgb(200, 200, 50)", #Yellow
    	#"rgb(230, 200, 10)", #Yellow
    	"rgb(255, 140, 0)", #orange?
    	"rgb(255, 140, 0)" #orange?
    ], inner=[2])
    t = (
               x=replace([0, 0, 1, 1, 0, 0, 1, 1], 0=> p1[1], 1=>p2[1]),
               y=replace([0, 1, 1, 0, 0, 1, 1, 0], 0=> p1[2], 1=>p2[2]),
               z=replace([0, 0, 0, 0, 1, 1, 1, 1], 0=> p1[3], 1=>p2[3]),
                i=[7, 0, 0, 0, 4, 4, 2, 6, 4, 0, 3, 7],
                j=[3, 4, 1, 2, 5, 6, 5, 5, 0, 1, 2, 2],
                k=[0, 7, 2, 3, 6, 7, 1, 2, 5, 5, 7, 6],
                facecolor=facecolor,
                flatshading=true

        )
    return t
end

function clear()
	deletetraces!(p,0:10)
end
