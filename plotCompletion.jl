using CSV, DataFrames, CairoMakie, StatsBase

# Generates X-coordinate array in a specific order (B,A,M,E).
# Returns dataframe with only data needed for plotting.
function calculateXColumn(df, order)
    xArr = Vector{Int}()
    grpByDomain = groupby(df, :dom)
    for grp in grpByDomain
        xPos = findfirst(==(grp.dom[1]), order)
        append!(xArr, fill(xPos, nrow(grp)))
    end
    domainDf = combine(grpByDomain, ["% completion", "% redundancy"])
    domainDf.x = xArr
    rename!(domainDf, ["% completion" => :comp, "% redundancy" => :red])
    return domainDf
end

function main(path)
    domainOrder = ["Bacteria", "Archaea", "Missing", "Eukarya"]
    colours = [Makie.wong_colors()[x] for x in [1,2,4,3]]
    binDB = coalesce.(CSV.read(path, DataFrame), "Missing")
    cm, pt = 96/2.54, 4/3
    plottingData = calculateXColumn(binDB, domainOrder)
    figure = Figure(fontsize = 14pt, size = (29.21cm, 12.09cm))
    # Subsections A and B.
    ax1 = Axis(figure[1,1], xlabel = "Genome completion (%)", ylabel = "Contamination (%)", title = "Anvi'o")
    ax2 = Axis(figure[1,2], xlabel = ax1.xlabel, ylabel = "SCG redundancy (%)", title = "EukCC2")
    for d in domainOrder
        tmpDomainDf = filter(:dom => ==(d), plottingData)
        idx = findfirst(==(d), domainOrder)
        d == "Eukarya" ? ax = ax2 : ax = ax1
        println(d, ": ", cor(tmpDomainDf.comp, tmpDomainDf.red))
        scatter!(ax, tmpDomainDf.comp, tmpDomainDf.red, alpha = 0.5, color = colours[idx])
    end
    linkxaxes!(ax1, ax2)
    # Subsection C.
    ax3 = Axis(figure[2, 1][1, 1:3], xticks = (1:3, domainOrder[1:3]), title = "Anvi'o", ylabel = "Genome completion (%)")
    ax4 = Axis(figure[2, 1][1, 4], xticks = ([4], ["Eukarya"]), title = "EukCC2")
    for d in domainOrder
        idx = findfirst(==(d), domainOrder)
        tmpBoxDf = filter(:dom => ==(d), plottingData)
        d == "Eukarya" ? ax = ax4 : ax = ax3
        boxplot!(ax, tmpBoxDf.x, tmpBoxDf.comp, color = colours[idx], strokewidth = 2)
    end
    linkyaxes!(ax3, ax4)
    hideydecorations!(ax4, grid = false)
    # Subsection D.
    domainFreqDf = combine(groupby(plottingData, :dom), :dom => countmap)
    domainFreqArr = [collect(values(x))[1] for x in domainFreqDf.dom_countmap]
    domainFreqDf.freqs = domainFreqArr
    println(domainFreqDf)
    ax5 = Axis(figure[2, 2], xticks = (1:4, domainOrder), ylabel = "Frequency")
    for x in eachrow(domainFreqDf)
        xPos = findfirst(==(x.dom), domainOrder)
        barplot!(ax5, xPos, x.freqs, strokewidth = 2, color = colours[xPos])
    end
    # Aligning Y-axis labels of A and C, and B and D to each other
    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax3])
    yspace2 = maximum(tight_yticklabel_spacing!, [ax5, ax2]) + 15
    ax1.yticklabelspace, ax3.yticklabelspace = Tuple(fill(yspace, 2))
    ax2.yticklabelspace, ax5.yticklabelspace = Tuple(fill(yspace2, 2))
    # Generates legend on the righthand side.
    Legend(figure[1:2, 3],
       [MarkerElement(marker = :circle, markercolor = colours[i], markersize = 20) for i in 1:4],
       domainOrder, "Domains", framevisible = false)
    # Generates subsection and N labels.
    subSecKwargs = (fontsize = 18pt, font = :bold, halign = :right, padding = (0, 5, 5, 0))
    nKwargs = (fontsize = 14pt, padding = (-100, 0, 0, 0), halign = :right)
    nA = nrow(filter(:dom => !=("Eukarya"), plottingData))
    nB = nrow(filter(:dom => ==("Eukarya"), plottingData))
    lblTbl = (x = vcat(fill(1,4), fill(2,4)), y = [1, 1, 2, 2, 1, 1, 2, 2],
          lbl = ["A", "n = $nA",
                 "B", "n = $nB",
                 "C","", "D", "n = $(nrow(plottingData))"])
    for i in 1:2:length(lblTbl.x)
        Label(figure[lblTbl.x[i], lblTbl.y[i], TopLeft()], lblTbl.lbl[i]; subSecKwargs...)
        Label(figure[lblTbl.x[i+1], lblTbl.y[i+1], TopRight()], lblTbl.lbl[i+1]; nKwargs...)
    end
    return figure
end

binPath = joinpath(@__DIR__, "private", "data", "outputDB.csv")
plot = main(binPath)
save(joinpath(@__DIR__, "private/", "plots", "completion.svg"), plot)
