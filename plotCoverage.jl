using CSV, DataFrames, CairoMakie

function main(bin, sample)
    binDB = coalesce.(CSV.read(bin, DataFrame), "Missing")
    sampleDB = CSV.read(sample, DataFrame)
    # Rounds number of samples (taken from sampleDB) to nearest multiple of 10. Used for histogram limits.
    roundedN = round(nrow(sampleDB) / 10) * 10
    samples = split.(binDB.sample, "-")
    sampleFreq = length.(samples)
    binLimit = round(maximum(sampleFreq) / 10) * 10
    # Plotting figure
    cm, pt = 96/2.54 , 4/3
    figure = Figure(fontsize = 14pt, size = (29.21cm, 12.09cm))
    # Subsection A.
    ax1 = Axis(figure[1, 1], xlabel = "Samples per bin", ylabel = "Frequency")
    hist!(ax1, sampleFreq, strokewidth = 2, bins=0:10:binLimit)
    # Subsection B.
    relAbs = [tryparse.(Float64, x) for x in split.(binDB.RA, "-")]
    ax2 = Axis(figure[1, 2], xlabel = "Samples per bin", ylabel = "Mean relative abundance (%)")
    scatter!(ax2,sampleFreq, mean.(relAbs))
    linkxaxes!(ax1,ax2)
    return figure
end

binPath = joinpath(@__DIR__, "private", "data", "outputDB.csv")
samplePath = joinpath(@__DIR__, "private", "data", "sampleData-v1.csv")
main(binPath, samplePath)