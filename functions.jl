# Various (somewhat) reusable functions.
using DataFrames, CSV, StringEncodings, CairoMakie, StatsBase
function constructDF(path)
    files = readdir(path, join = true)
    # Below line makes sure we're not trying to read files like .DS_Store
    nonDotFiles = filter(f -> !startswith(basename(f), "."), files)
    nFiles = length(nonDotFiles)
    # Pre-allocate arrays
    bins = Vector{Vector{String}}(undef, nFiles)
    relAbundances = Vector{Vector{Float64}}(undef, nFiles)
    samples = Vector{String}(undef, nFiles)
    # Below FOR: record coverage data of for all bins >0 relataive abundance for all bins
    # This is stored in nested arrays - each nested array represents one file (sample)
    for (x, file) in enumerate(nonDotFiles)
        tmpDF = CSV.read(file, DataFrame)
        filter!(col -> col[2] > 0, tmpDF)
        bins[x] = tmpDF[!, 1]
        relAbundances[x] = tmpDF[!, 2]
        samples[x] = basename(file)
    end
    # Pre-allocate sample array for FOR loop
    finalS = Vector(undef, sum(length.(bins)))
    endIdx = 1
    # Expand sample array to same size as flattened bins/RAs arrays
    # If sample1 has 300 bins in it, we need finalS[1:300] = "sample1"
    for idx in eachindex(samples)
        nPerSample = length(bins[idx])
        startIdx = endIdx
        endIdx += nPerSample - 1
        finalS[startIdx:endIdx] .= samples[idx]
        endIdx += 1
    end
    return DataFrame(bin = reduce(vcat, bins), sample = finalS, RA = reduce(vcat, relAbundances))
end

# Filtering rules for bins provided by Lisa.
function filterBins(col)
    notBinStart = .!startswith.(col, "bin.")
    noMetabat = .!contains.(col, "metabat")
    noMaxBin = .!contains.(col, "maxbin")
    mappedOnly = .!=(col, "unmapped")
    notBinStart && noMetabat && noMaxBin && mappedOnly
end

# Remove unecessary sections from PID names.
function formatPIDs(arr)
    fmtArr = []
    splitArr = split.(arr, "_")
    # Below pop!() removes the "_result" from end of PID.
    pop!.(splitArr)
    # Iteratively removes channel number where necessary from all PIDs.
    # Runs fmtLastElem() before pushing corrected name to array.
    for arr in splitArr
        if contains(arr[end], "S")
            out = join(arr[1:end-1], "_")
            push!(fmtArr, rmLChr(out))
        else
            out = join(arr, "_")
            push!(fmtArr, rmLChr(out))
        end
    end
    return fmtArr
end

# Checks to if last PID chr is a letter, if so remove it.
rmLChr(pid) = isletter(pid[end]) ? pid[1:end-1] : pid

# Remove all run 2 (PID_1100_...) instances of duplicated PIDs in samples.
function rmDuplicatePIDs!(df, col)
    uniquePIDs = unique(df[!, col])
    splitPIDs = split.(uniquePIDs, "_")
    sampleIDs = getindex.(splitPIDs, 3)
    idFreqs = countmap(sampleIDs)
    duplicatedIDs = collect(keys(filter(x -> last(x) > 1, idFreqs)))
    duplicatedPIDs = "PID_1100_" .* duplicatedIDs
    filter!(:sample => s -> !(s in duplicatedPIDs), df)
end

rmRunPrefix(arr) = [join(x[2:end], "_") for x in split.(arr, "_")]

function rmDuplicateBins!(df, col, path)
    mergedData = CSV.read(path, DataFrame)
    mergedBins = collect(Iterators.flatten(split.(mergedData.merged, ",")))
    # mergedBins contains representatives (ones to keep), need to find the other ones.
    binsToDrop = setdiff(mergedBins, mergedData.representative)
    # One-liner to remove "runN" from bins.
    fmtdBTD = rmRunPrefix(binsToDrop)
    filter!(col => !(in(fmtdBTD)), df)
end

function transformDF(df)
    nRows = length(unique(df.bin))
    binArr = Vector{String}(undef, nRows)
    samplesArr = Vector{String}(undef, nRows)
    RAsArr = Vector{String}(undef, nRows)
    grpByBins = groupby(df, :bin)
    for (i, grp) in enumerate(grpByBins)
        binArr[i] = grp.bin[1]
        samplesArr[i] = join(grp.sample, "-")
        RAsArr[i] = join(grp.RA, "-") 
    end
    return DataFrame(bin = binArr, sample = samplesArr, RA = RAsArr)
end

function findNonEukPhyla(df, dom)
    tmpDf = filter(:dom => ==(dom), df)
    taxa = split.(tmpDf.tax, "_")
    phyla = countmap(getindex.(taxa, 2))
    outDf = DataFrame(phyla = collect(keys(phyla)), freq = collect(values(phyla)))
    sort!(outDf, :phyla)
end