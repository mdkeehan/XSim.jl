function getOurPhenVals(my::Cohort, varRes=0)
    if varRes == 0
       if length(common.LRes) != 0
           LRes = common.LRes
        else
           error("getOurPhenVals(): common.LRes is not initialized")
        end
    else
       if typeof(varRes)==Float64
          LRes = cholesky(fill(varRes,1,1)).U'
       else
          LRes = cholesky(varRes).U'
        end
    end

    n = size(my.animalCohort,1)
    nTraits = size(common.LRes,2)
    phenVals = Array{Float64,2}(undef,n,nTraits)
    genVals  = getOurGenVals(my,nTraits)
    for (i,animal) = enumerate(my.animalCohort)
        if length(animal.phenVal) == 0
	    #for t in 1:nTraits?
            animal.phenVal =  animal.genVal + LRes * randn(nTraits)
        end
        phenVals[i,:] = animal.phenVal
    end
    return phenVals
end

function getOurGenVals(my::Cohort,nTraits=0)
    if nTraits == 0
        if  size(common.LRes,2) != 0           
            nTraits = size(common.LRes,2)
        else
            error("getOurGenVals(): common.LRes is not initialized, cannot get the number of traits from the 2nd dimension of the residual variance matrix.")
        end
    else
	print("nTraits = $(nTraits)\n")
	print("qtl_effects = $(common.G.qtl_effects)\n")
	print("qtl_index = $(common.G.qtl_index)\n")
    end
    #if nTraits != 1
	#error("nTraits == $(nTraits) is not yet implmented")
    #end
     
    n = size(my.animalCohort,1)
    genVals = Array{Float64}(undef,n,nTraits)
    traitIndex = 1 # hardwire only 1 trait at the moment
    for (i,animal) = enumerate(my.animalCohort)
        if i%1000 == 0
            println("getOurGenVals(): ", i)
        end
        if length(animal.genVal) == 0
            getMyHaps(animal)
            myGenotypes = getMyGenotype(animal)
            #animal.genVal = dot(myGenotypes[common.G.qtl_index],common.G.qtl_effects)
	    # need something like
            #for traitIndex in 1:nTraits ?
#print("myGenotypes[common.G.qtl_index] : type $(typeof(myGenotypes[common.G.qtl_index])) value $(myGenotypes[common.G.qtl_index])\n")
#print("animal.genVal : type $(typeof(animal.genVal)) value $(animal.genVal)\n")
              push!(animal.genVal, dot(myGenotypes[common.G.qtl_index],common.G.qtl_effects))
	    #end
        end
#	print("animal.genVal = $(animal.genVal)\n")
        genVals[i,traitIndex] = animal.genVal[traitIndex]
    end
    return genVals
end
