using CSV, CairoMakie, DataFrames

mx = Matrix(CSV.read(joinpath(@__DIR__, "private", "data", "concurrencyMatrix.csv"), DataFrame))
#heatmap(mx)
binDB = coalesce.(CSV.read(joinpath(@__DIR__, "private/data/outputDB.csv"), DataFrame), "Missing")
mxIndex = CSV.read(joinpath(@__DIR__, "private", "data", "concurrencyVector.csv"), DataFrame)[!, 1]

domIdxs = Dict()
grpByDom = groupby(binDB[:, [1,end]], :dom)
# Find indexes of all memebrs of each domain
[domIdxs[grp[1,2]] = indexin(grp.bin, mxIndex) for grp in grpByDom]

dataPerDom = Dict()
domains = "Bacteria", "Archaea", "Missing", "Eukarya"



#[dataPerDom[dom] = mx[:, domIdxs[dom]] for dom in domains]
plotTbl = (x = [], y = [])
for (i, dom) in enumerate(domains)
    tmpData = filter(x -> x > 0.5, vec(mx[:, domIdxs[dom]]))
    append!(plotTbl.x, fill(i, length(tmpData)))
    append!(plotTbl.y, tmpData)
end
barplot(plotTbl.x, plotTbl.y)