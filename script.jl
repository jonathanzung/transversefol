function setup(isosig)
	include("batch/$(isosig).txt")
	global bt=BoundaryTriangulation(fans, face_coorientations, firstrungs,alledges,weights)
	#=
	global longitude
	if longitude == nothing
		longitude=find_longitude(fans)#weights of the different faces
	end
	=#
	global longitudes = find_longitudes_iterative(fans,1000)
	#push!(longitudes, longitude)
	#sort!(longitudes, by=l->all_slopes(longitude_to_candidate(bt,l),time=30000)[1])
	
	global p=PlotlyJS.plot()
	#=
	function valid_slope(s)
		return true
		#return maximum(abs.(s)) <= 5 && s[1] <= 0
	end
	filter!(l-> valid_slope(all_slopes(longitude_to_candidate(bt,l),time=2000)), longitudes)
	=#

	global Elong=Envelope()
	global Econstr=PEnvelope()
	for l in longitudes
		local c
		c=longitude_to_candidate(bt,l)
		x=approximant_all_slopes(c, time=400*sum(l))
		push!(Elong, (x,c))
		for constr in constraints(Longitude(bt,l))
			push!(Econstr, (constr, c))
		end
	end
	global upperE = Envelope{Upper}(copy(Elong.A))
	global lowerE = Envelope{Lower}(copy(Elong.A))
	add_trace!(p, _plotjs(Elong))
	add_trace!(p, _plotjs(Econstr))

	display(p)

	println("done setup")
end

include("find_surface.jl")
include("batch/2cusp_manifest.txt")
isosig = "siddhi2"
#isosig = "eLMkbcddddedde_2100"
#isosig = "gvLQQcdeffeffffaafa_201102"
#isosig = "gLLAQcdecfffhsermws_122201"
#isosig = "fLLQcbecdeepuwsua_20102"
#isosig = "fLLQcbeddeehhbghh_01110"
#isosig = "challenge2"
#isosig = "gLLPQbefefefhhhhhha_011102"
#isosig = "gLLPQcdfefefuoaaauo_022110"
#isosig = "gLLPQcdfefefuoaaauo_022110"
#setup(isosig)




function regimen(E::Envelope{S}) where {S}
	if S==Upper
		#target=[-1.35,0.075]
		target=[10,10]
	elseif S==Lower
		target=[-10,-10]
	else
		@assert false
	end
	E = try_improve(E; nsubdivide=1, iters=30000, time=1000, target=target, radius=0.002)
	add_trace!(p, _plotjs(E))
	E = try_improve(E; nsubdivide=0, iters=100000, time=2000, target=target, radius=0.001, beta=800)
	add_trace!(p, _plotjs(E))
	E = try_improve(E; nsubdivide=0, iters=1000000, time=2000, target=target, radius=0.001, beta=1500)
	add_trace!(p, _plotjs(E))
	return E
end




#=
dummy_candidate=random_candidate(bt,0)
L6a2E=Envelope()
push!(L6a2E, (T[0,0], dummy_candidate))
for pt in [(-4,1/2+0.000001), (-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
	push!(L6a2E, (T[pt...], dummy_candidate))
end
=#



randE = random_trials(bt)
add_trace!(p, _plotjs(randE))
#display(p)
#add_trace!(p, _plotjs(E2, maxabs=40))
#

using ProfileView
upperE = regimen(upperE)
Profile.init(n=2*10^7, delay=0.005)
lowerE = @profile regimen(lowerE)
ProfileView.view()

addtraces!(p, _plotjs(lowerE, upperE)...)


#add_trace!(p, _plotjs(L6a2E, fill=true))
#add_trace!(p, _plotjs(accurate_E,fill=true))
#scatter_envelope!(p, accurate_E)
#scatter!(p, [-2,-1,-1/2,0],[1/2,1/3,1/6,0])
#scatter!(p,[-2,-1,-1/2],[1/3,1/6,0])


#=
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
=#


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

#(0,0) is triangulated with 4 ideal triangles, has at least 2 punctures
# V - 3/2 *4 + 4 = V - 6 + 4 = V - 2. So V=2 => torus, V=4 => sphere. It's either a twice punctured torus or a 4 times punctured sphere 
#
#
#



# m125
# Records
# [[0.5, -2.0], [0.014492753623188406, -1.3714285714285714], [-1.375, 0.06666666666666667], [0.25, -1.5], [-1.3703703703703705, 0.014492753623188406], [0.0136986301369863, -1.368421052631579], [-2.0, 0.5], [-1.4, 0.1111111111111111], [0.1111111111111111, -1.4], [0.07142857142857142, -1.375], [0.0, 0.0], [-1.5, 0.25]]

