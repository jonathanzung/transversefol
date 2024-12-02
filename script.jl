using Serialization
using DataFrames

include("search.jl")
include("find_surface.jl")

function setup(isosig; nlongs = 30)
	println("setting up $(isosig)")
	flush(stdout)
	include("batch/$(isosig).txt")


	global bt=BoundaryTriangulation(fans, face_coorientations, firstrungs,alledges,rungs)

	#=
	global longitude
	if longitude == nothing
		longitude=find_longitude(fans)#weights of the different faces
	end
	=#
	#global longitudes = find_longitudes_iterative(fans,1000)
	#push!(longitudes, longitude)
	#sort!(longitudes, by=l->all_slopes(longitude_to_candidate(bt,l),time=30000)[1])
	
	global ncusps = length(bt.firstrungs)
	
	#=
	function valid_slope(s)
		return true
		#return maximum(abs.(s)) <= 5 && s[1] <= 0
	end
	filter!(l-> valid_slope(all_slopes(longitude_to_candidate(bt,l),time=2000)), longitudes)
	=#

	global Elong=PEnvelope()

	global longitudeDF = DataFrame()
	global long_dict = DefaultDict(()->[])

	for l in find_longitudes_iterative(fans)
		#ss=approximant_all_slopes(c, time=400*sum(l))
		ss = [y//x for (x,y) in slopes(Longitude(bt,l))]

		if !haskey(long_dict, ss)
			@show length(long_dict)
			flush(stdout)
		end
		push!(long_dict[ss], l)


		if length(long_dict) >= nlongs
			break
		end
	end

	for (ss, ls) in long_dict
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		push!(Elong, (ss, c))
		for (i,constr) in enumerate(constraints(Longitude(bt,ls[i])))
			push!(Econstr[i], (constr, c))
		end
		push!(longitudeDF, (ss=slopes(Longitude(bt,ls[i])), l=ls[i], chi = -sum(ls[i])//2, normalizedchi = normalizedchi(Longitude(bt,ls[i]))))
	end

	global long_ranges = [[minimum(filter(r -> abs(r) < CLIP, collect(x[i] for x in keys(long_dict)))),
						   maximum(filter(r -> abs(r) < CLIP, collect(x[i] for x in keys(long_dict))))] for i in 1:ncusps]
	global p=PlotlyJS.plot()

	global Eupper = Envelope{Upper}(copy(Elong.A))
	global Elower = Envelope{Lower}(copy(Elong.A))


	println("done setup")
	flush(stdout)
end
function regimen(E::Envelope{S}) where {S}
	if S==Upper
		target=[CLIP,CLIP]
	elseif S==Lower
		target=[-CLIP,-CLIP]
	else
		@assert false
	end
	E = try_improve(E; nsubdivide=1, iters=30000, time=1000, target=target, radius=0.002)
	flush(stdout)
	if isinteractive()
		add_trace!(p, _plotjs(E))
	end
	E = try_improve(E; nsubdivide=0, iters=100000, time=2000, target=target, radius=0.001, beta=800)
	flush(stdout)
	if isinteractive()
		add_trace!(p, _plotjs(E))
	end
	E = try_improve(E; nsubdivide=0, iters=1000000, time=2000, target=target, radius=0.001, beta=1300)
	flush(stdout)
	if isinteractive()
		add_trace!(p, _plotjs(E))
	end
	return E
end

function runjob(i)
	include("batch/2cusp_manifest.txt")
	setup(isosig)
	global p
	#isosig=isosigs[i]
	#setup(isosig)

	#randE = random_trials(bt)
	#add_trace!(p, _plotjs(randE))

	global Eupper = regimen(Eupper)
	global Elower = regimen(Elower)
	addtraces!(p, _plotjs(Elower, Eupper)...)
	add_trace!(p, _plotjs(Elong; color=LONGITUDE_COLOUR))
	for i in 1:4
		add_trace!(p, _plotjs(Econstr[i]))
	end

	xe = 0.5 * (long_ranges[1][2]-long_ranges[1][1])
	ye = 0.5 * (long_ranges[2][2]-long_ranges[2][1])
	update_xaxes!(p,range=[long_ranges[1][1]-xe, long_ranges[1][2]+xe],autorange=false, title="cusp 1 slope")
	update_yaxes!(p,range=[long_ranges[2][1]-ye, long_ranges[2][2]+ye], autorange=false, title="cusp 2 slope")

	PlotlyJS.savefig(p, "batch/$(isosig).html")
	serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong, longitudeDF=longitudeDF))
	flush(stdout)
end

function quickview(i)
	include("batch/2cusp_manifest.txt")
	#isosig=isosigs[i]
	isosig="ivLLQQccfhfeghghwadiwadrv_20110220"
	#isosig="hLLLQkbeegefgghhhahabg_0111022"
	isosig="jLvvQQQbhigghihgixaxxvvvvcc_102222010"
	setup(isosig)




	global longitudeDF = DataFrame()
	for (ss, ls) in long_dict
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		push!(longitudeDF, (ss=slopes(Longitude(bt,ls[i])), l=ls[i], chi = -sum(ls[i])//2, normalizedchi = normalizedchi(Longitude(bt,ls[i]))))
	end

	global long_ranges = [[minimum(filter(r -> abs(r) < CLIP, collect(x[i] for x in keys(long_dict)))),
						   maximum(filter(r -> abs(r) < CLIP, collect(x[i] for x in keys(long_dict))))] for i in 1:ncusps]

	config = PlotConfig(modeBarButtonsToAdd=[
    "drawline",
    "drawopenpath",
    "drawclosedpath",
    "drawcircle",
    "drawrect",
    "eraseshape"
	])
	global p=PlotlyJS.plot(config=config)

	global Eupper = Envelope{Upper}(copy(Elong.A))
	global Elower = Envelope{Lower}(copy(Elong.A))

	Elower2, Eupper2 = extreme_candidates(bt)



	#global p=PlotlyJS.plot()
	global Econstr=[PEnvelope() for i in 1:4]

	#display(p)
	#isosig=isosigs[i]
	#setup(isosig)

	#randE = random_trials(bt)
	#add_trace!(p, _plotjs(randE))

	#global Eupper = regimen(Eupper)
	#global Elower = regimen(Elower)
	#global tup = deserialize("/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls")
	addtraces!(p, _plotjs(Elower, Eupper)...)

	addtraces!(p, _plotjs(Elower2, Eupper2)...)
	#add_trace!(p, _plotjs(Elong; color=LONGITUDE_COLOUR))

	global longitudeDF = DataFrame()
	for (ss, ls) in long_dict
		_, i = findmin(x-> (count(y->y==0, x), sum(x.^2)), ls)
		c=longitude_to_candidate(bt,ls[i])
		nchi = normalizedchi(Longitude(bt,ls[i]))
		push!(longitudeDF, (ss=slopes(Longitude(bt,ls[i])), x=ss[1], y=ss[2], l=ls[i], chi = -sum(ls[i])//2, nchi = nchi, text=string((nch=nchi,ss=ss))))
		for (i,constr) in enumerate(constraints(Longitude(bt,ls[i])))
			push!(Econstr[i], (constr, c))
		end
	end
	add_trace!(p, PlotlyJS.scatter(longitudeDF,x=:x, y=:y, marker=attr(line=attr(width=0), size=25 ./ log.(4 .- longitudeDF[!,:nchi]), color=LONGITUDE_COLOUR), text=:text, mode="markers"))
	for i in 1:4
		add_trace!(p, _plotjs(Econstr[i]))
	end


	xe = 0.5 * (long_ranges[1][2]-long_ranges[1][1])
	ye = 0.5 * (long_ranges[2][2]-long_ranges[2][1])
	update_xaxes!(p,range=[long_ranges[1][1]-xe, long_ranges[1][2]+xe],autorange=false, title="cusp 1 slope")
	update_yaxes!(p,range=[long_ranges[2][1]-ye, long_ranges[2][2]+ye], autorange=false, title="cusp 2 slope")



	#PlotlyJS.savefig(p, "batch/$(isosig).html")
	#serialize("batch/$(isosig).jls", (bt=bt, Eupper=Eupper, Elower=Elower, Elong=Elong))
	flush(stdout)
	p
end

function aggregate_bounds()
	include("batch/2cusp_manifest.txt")
	#isosig=isosigs[i]
	isosig="ivLLQQccfhfeghghwadiwadrv_20110220"
	#isosig="hLLLQkbeegefgghhhahabg_0111022"
	#setup(isosig)
	#
	


	for isosig in isosigs[30:50]
	end

	global p=PlotlyJS.plot(config=config)

	global Eupper = Envelope{Upper}(copy(Elong.A))
	global Elower = Envelope{Lower}(copy(Elong.A))

	Elower2, Eupper2 = extreme_candidates(bt)



	#global p=PlotlyJS.plot()
	global Econstr=[PEnvelope() for i in 1:4]

	#display(p)
	#isosig=isosigs[i]
	#setup(isosig)

	#randE = random_trials(bt)
	#add_trace!(p, _plotjs(randE))

	#global Eupper = regimen(Eupper)
	#global Elower = regimen(Elower)
	global tup = deserialize("/home/jonathan/engaging_sshfs/transversefol/batch/$(isosig).jls")
end

#isosig = "siddhi2"
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

#=
isosig = "eLMkbcddddedde_2100"
setup(isosig)

Profile.init(n=10^7, delay=0.01)
@profile setup(isosig)
using ProfileView
if isinteractive()
	ProfileView.view()
end
=#



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




#=
dummy_candidate=random_candidate(bt,0)
L6a2E=Envelope()
push!(L6a2E, (T[0,0], dummy_candidate))
for pt in [(-4,1/2+0.000001), (-2,1/2), (-1, 1/3), (-1/2, 1/6), (-1/3, 1/9), (-1/4, 1/12), (-1/5, 1/15), (-1/6, 1/18)]
	push!(L6a2E, (T[pt...], dummy_candidate))
end
=#
