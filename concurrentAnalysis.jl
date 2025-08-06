using CSV, StatsBase, DataFrames, Tables
binDB = coalesce.(CSV.read(joinpath(@__DIR__, "private/data/outputDB.csv"), DataFrames.DataFrame), "Missing")
binDB.sample = split.(binDB.sample, "-")
samples = unique(collect(Iterators.flatten(binDB.sample)))
bins = unique(binDB.bin)

function generateSampleArray(sampleArr, df)
    out = Vector(undef, length(sampleArr))
    for (i, s) in enumerate(sampleArr)
        tmp = filter(:sample => x -> s in x, df)
        out[i] = tmp.bin
    end
    return out
end

function generateMx(uniqueBins, binVec)
    n = length(uniqueBins)
    concurrencyMatrix = Matrix{Float64}(zeros(n, n))
    maxArr = Vector(undef, n)
    for (k, z) in enumerate(uniqueBins)
        tmpArrs = filter(x -> z in x, binVec)
        associatedBins = collect(Iterators.flatten(tmpArrs))
        maxArr[k] = length(tmpArrs)
        for y in associatedBins
            i = findall(q -> q == z, uniqueBins)[1]
            j = findall(q -> q == y, uniqueBins)[1]
            concurrencyMatrix[i, j] += 1
        end
    end
    return concurrencyMatrix ./ maxArr
end

binsPerSample = generateSampleArray(samples, binDB)
mx = generateMx(bins, binsPerSample)
outPath = joinpath(@__DIR__, "private","data","concurrencyMatrix.csv")
CSV.write(outPath, Tables.table(mx))