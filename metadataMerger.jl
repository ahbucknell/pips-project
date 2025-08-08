# Merges various metadata files togethe into some sample database.
using CSV, DataFrames
include("functions.jl")
function main()
    binPath = joinpath(@__DIR__, "private", "data", "outputDB.csv")
    binDB = CSV.read(binPath, DataFrame)
    allSamples = split.(binDB.sample, "-")
    samples = unique(collect(Base.Iterators.flatten(allSamples)))
    sampleDB = DataFrame(sample = samples)
    # Merging watersheds to PIDs from binDB.
    watershedPath = joinpath(@__DIR__, "private", "data", "prerequisite", "watersheds.txt")
    watersheds = CSV.read(watershedPath, DataFrame, header = false)
    ## Removes "+j" section of "i+j" suffix found on some PIDs, where i and j are letters.
    rmSuffix = getindex.(split.(watersheds[!, 1], "+"), 1)
    ## Removes the remaining "i" section of the "i+j" suffix - PIDs are now clean.
    watersheds[!,1] = rmLChr.(rmSuffix)
    rename!(watersheds, ["sample", "watershed"])
    sampleDF = outerjoin(sampleDB, watersheds, on = :sample)

    # Merging misc. metadata
    metadataPath = joinpath(@__DIR__, "private", "data", "prerequisite", "sampleMetadata.txt")
    metadataDF = CSV.read(metadataPath, DataFrame)
    metadataDF.sample = rmLChr.(metadataDF.sample)
    ## Remove metadata rows where sample is "N" or "NA".
    filter!(:sample => n -> !(occursin("N", n)), metadataDF)
    ## Merge all metadataDF rows to sampleDF rows - does result in some missings.
    sampleDF = leftjoin(sampleDF, metadataDF, on = :sample)

    # Merging sample types
    sampletypePath = joinpath(@__DIR__, "private", "data", "prerequisite", "sampleTypes.txt")
    sampletypeDF = CSV.read(sampletypePath, DataFrame)
    rename!(sampletypeDF, :Sequenced_sample => :sample)
    select!(sampletypeDF, ["sample", "SampleType"])
    sampleDF = leftjoin(sampleDF, sampletypeDF, on = :sample)

    outPath = joinpath(@__DIR__, "private/", "data/", "sampleData-base.csv")
    CSV.write(outPath, sampleDF)
end
main()