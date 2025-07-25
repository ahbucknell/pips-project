# Merges various metadata files together.
using CSV, DataFrames
include("functions.jl")


binPath = joinpath(@__DIR__, "private","data","outputDB.csv")
binDB = CSV.read(binPath, DataFrame)
allSamples = split.(binDB.sample, "-")
samples = unique(collect(Base.Iterators.flatten(allSamples)))
sampleDB = DataFrame(sample = samples)
# Merging watersheds to PIDs from binDB.
watershedPath = joinpath(@__DIR__, "private","data", "prerequisite", "watersheds.txt")
watersheds = CSV.read(watershedPath, DataFrame, header = false)
## Removes "+j" section of "i+j" suffix found on some PIDs, where i and j are letters.
rmSuffix = getindex.(split.(watersheds[!, 1], "+"), 1)
## Removes the remaining "i" section of the "i+j" suffix - PIDs are now clean.
watersheds[!,1] = rmLChr.(rmSuffix)
rename!(watersheds, ["sample", "watershed"])
sampleDF = innerjoin(sampleDB, watersheds, on = :sample)




