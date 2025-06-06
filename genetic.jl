function genetic(f, population::Vector{T}, jiggle, crossover, nsteps) where {T}
    newpopulation = Vector{T}
 
    @threads for c in sample(population, nmutate)
        push!(newpopulation, jiggle(c))
    end
    @threads for c in sample(population, nbreed)
        c2 = rand(population)
        push!(newpopulation, crossover(c,c2))
    end
    population = cull(newpopulation, fitness, popsize)
end

function cull(population, fitness, popsize)
end

function crossover(c1::Cand{DiscreteHomeo}, c2::Cand{DiscreteHomeo}; prob=0.2)
    @assert c1.bt == c2.bt
    struct Cand{S<:Homeo} #should it be mutable?
        bt::BoundaryTriangulation
        d::ArrayDict{Junction, S, 2}
    end

    d3=copy(c1.d)
    for i in eachindex(d3)
        d3[i] = c2.d[i]
    end

    return Cand{DiscreteHomeo}(bt, d3)
end



function annealing(f, initial, jiggle, betastart, betafinish, nsteps; verbose=true, minacc = 100, maxacc = 5000)
	#linear annealing on range betastart, betafinish
	current = initial
	curracc = minacc
	currval = f(current; acc=curracc)
	reject_count = 0
	accept_count = 0

	vals=[]
	push!(vals,currval)
	for (i,currbeta) in zip(1:nsteps, range(betastart,betafinish, nsteps))
		jig = jiggle(current)
		newacc = minacc
		newval = f(jig; acc=minacc)

		r=rand()

		@label here
		dE = exp(currbeta * (-newval+ currval))
		prob = 1/(1+dE)
		#prob will be an interval
	
		#@show prob
		if r > prob
			reject_count += 1
		elseif r < prob
			current = jig
			currval = newval
			curracc = newacc
			accept_count += 1
		elseif (curracc >= maxacc && newacc >= maxacc)
			println("not enough accuracy")
			@show prob, dE
			@show currval
			@show newval
	#		@assert false
		else #improve the accuracy of our computation
			if newacc < maxacc
				newacc *= 3
				newval = f(jig; acc=newacc)
			end
			if curracc < newacc
				curracc *= 3
				currval = f(current; acc = curracc)
			end
			@goto here
		end
		if verbose
			push!(vals, currval)
		end

		if verbose && i%10000 == 0
			@show (exact_slope(current))
			@show (accept_count, reject_count) 
			p=plot(1:length(vals), vals)
			display(p)
		end
	end
	@show (accept_count, reject_count)
	return current
end
