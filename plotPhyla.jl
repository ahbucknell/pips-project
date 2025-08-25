using CSV, DataFrames, CairoMakie, StatsBase, Taxonomy
include("functions.jl")

function findEukaryoticPhyla(df)
    eukDB = filter(:dom => ==("Eukarya"), df)
    fullTaxonomy = split.(eukDB.tax, "_")
    allLineages = Vector(undef, length(fullTaxonomy))
    # Convert taxonomy IDs to full taxonomy lineages.
    for (i, taxon) in enumerate(fullTaxonomy)
        finalTaxonID = tryparse(Int, taxon[end])
        lineage = Lineage(Taxon(finalTaxonID))
        allLineages[i] = lineage
    end
    phyla = Vector(undef, length(allLineages))
    # Not all the lineages have a phylum - if they don't and they belong to SAR, find the taxon after SAR.
    # If they don't belong to SAR and don't have a phyla, then count as Missing if the lineage length < 4, otherwise take the 5th element.
    for (j ,lineage) in enumerate(allLineages)
        try
            phyla[j] = name(lineage[:phylum])
        catch
            if "Sar" in name.(lineage)
                sarIdx = findfirst(==("Sar"), name.(lineage))
                phyla[j] = name(lineage[sarIdx + 2])
            else
                length(lineage) < 4 ? phyla[j] = "Missing" : phyla[j] = name(lineage[4])
            end
        end
    end
    return phyla
end

function main(path)
    domainOrder = ["Bacteria", "Archaea", "Eukarya"]
    # If not already loaded, load the Taxonomy.jl database.
    taxonomyPrivatePath = joinpath(@__DIR__, "private", "taxonomyPath.txt")
    taxonomyDBpath = readlines(taxonomyPrivatePath)[1]
    try
        current_db()
    catch
        taxDB = Taxonomy.DB(taxonomyDBpath * "/nodes.dmp", taxonomyDBpath * "/names.dmp")
    end
    binDB = coalesce.(CSV.read(path, DataFrame), "Missing")
    # Generate DF of eukaryotic phyla, move Missing to end of DF.
    eukPhylaDict = countmap(findEukaryoticPhyla(binDB))
    eukPhylaDf = DataFrame(phyla = collect(keys(eukPhylaDict)), freq = collect(values(eukPhylaDict)))
    println(eukPhylaDf)
    sort!(eukPhylaDf, :phyla)
    ePhylaDf = filter(:phyla => !=("Missing"), eukPhylaDf)
    missingOnly = filter(:phyla => ==("Missing"), eukPhylaDf)
    append!(ePhylaDf, missingOnly)
    # Generate non-eukaryotic phyla frequency DFs.
    bPhylaDf = findNonEukPhyla(binDB, "Bacteria")
    aPhylaDf = findNonEukPhyla(binDB, "Archaea")
    # Plotting figure
    colours = [Makie.wong_colors()[x] for x in [1,2,3]]
    cm, pt = 96/2.54, 4/3
    figure = Figure(fontsize = 14, size = (29.21cm, 12.09cm))
    dataOrder = [bPhylaDf, aPhylaDf, ePhylaDf]
    for (x, i) in enumerate([(1,1:3), (2,3), (2,1:2)])
        ax = Axis(figure[first(i), last(i)], xticks = (1:nrow(dataOrder[x]), dataOrder[x].phyla),
                  xticklabelrotation=pi/6, yscale = log10, yticks = LogTicks(0:3), ylabel = "Frequency",
                  title = domainOrder[x])
        barplot!(ax, 1:nrow(dataOrder[x]), dataOrder[x].freq, strokewidth = 1, color = colours[x])
        ylims!(0.5, 1000)
    end
    return figure
end

binPath = joinpath(@__DIR__, "private", "data", "outputDB.csv")
plot = main(binPath)
save(joinpath(@__DIR__, "private/", "plots", "phyla.svg"), plot)