using CSV, DataFrames, CairoMakie, StatsBase
include("functions.jl")

function main()
    binPath = joinpath(@__DIR__, "private", "data", "outputDB.csv")
    binDB = coalesce.(CSV.read(binPath, DataFrame), "Missing")
    highQualDB = filter(row -> row."% completion" >= 95 && row."% redundancy" <= 5, binDB)

    nPhylaDF = DataFrame(phyla = [], freq = [], dom = [])
    domainOrder = ["Bacteria", "Archaea", "Missing", "Eukarya"]
    for dom in domainOrder
        tmpDf = findNonEukPhyla(highQualDB, dom)
        tmpDf.dom .= dom
        nPhylaDF = vcat(nPhylaDF, tmpDf)
    end
    sort!(nPhylaDF, [:dom, :phyla], rev = [true, false])

    figure = Figure()
    
    colours = [Makie.wong_colors()[x] for x in [1,2,4,3]]
    colourDict = Dict(domainOrder .=> colours)

    n1 = nrow(nPhylaDF)
    ax1 = Axis(figure[1,1], xticks = (1:n1, nPhylaDF.phyla), xticklabelrotation = 45, ylabel = "Bint count")
    barplot!(ax1, 1:n1, nPhylaDF.freq, strokewidth = 2, color = [colourDict[dom] for dom in nPhylaDF.dom])

    highQualDB
end
main()