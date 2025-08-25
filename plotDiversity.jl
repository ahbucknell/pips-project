using CairoMakie, StatsBase, DataFrames, CSV
binPath = joinpath(@__DIR__, "private","data","outputDB.csv")
samplePath = joinpath(@__DIR__, "private","data","sampleData-v1.csv")
binDB = coalesce.(CSV.read(binPath, DataFrame), "Missing")
sampleDB = CSV.read(samplePath, DataFrame)
nonMissing = filter(:dom => !=("Missing"), binDB)

# Find the first instance of "None" and remove all elements after it per array.
function reFmtTaxArray(bigArr)
  noNones = Vector(undef, length(bigArr))
  for (x, arr) in enumerate(bigArr)
    if in("None", arr)
      noneIdx = findfirst(x -> x == "None", arr)
      tmpArr = arr[1:noneIdx-1]
      noNones[x] = tmpArr
    else
      noNones[x] = arr
    end
  end
  return noNones
end

function generateShannonIndex(arr)
  generaFreq = countmap(arr)
  freqs = collect(values(generaFreq))
  total = sum(freqs)
  props = freqs ./ total
  sdi = -(sum(props .* log.(props)))
  return sdi
end

# Converts a dictionary of watershed data into table for plotting
function dictToPlotArr(dict)
  tmpData = (x = [], y = [], z = [])
  for (x,(k,v)) in enumerate(dict)
    append!(tmpData.x, fill(x, length(v)))
    append!(tmpData.y, v)
    push!(tmpData.z, k)
  end
  return tmpData
end

function main()
    taxonomy = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    anvioData = filter(:dom => !=("Eukarya"), nonMissing)
    splitProkTaxa = split.(anvioData.tax, "_")
    refmtdProkTax = reFmtTaxArray(splitProkTaxa)
    anvioData[!, :fmtdTax] = refmtdProkTax
    atLeastGenus = filter(:fmtdTax => n -> length(n) >= 6,anvioData)
    atLeastGenus.s = split.(atLeastGenus.sample, "-")
    atLeastGenus[!, :genus] = getindex.(atLeastGenus.fmtdTax, 6)
    allSDIs = Dict{String, Float64}()
    # Un-nest the sample array and find only unique samples
    samples = unique(collect(Iterators.flatten(atLeastGenus.s)))
    for sample in samples
        perSampleDF = filter(:s => s -> in(sample, s), atLeastGenus)
        sampleSDI = generateShannonIndex(perSampleDF.genus)
        allSDIs[sample] = sampleSDI
    end
    wshedElevs, wshedSDIs = Dict(), Dict()
    plottingData =  DataFrame()
    sampleWsheds = groupby(sampleDB, :wshed)
    for grp in sampleWsheds
        currWshed = grp.wshed[1]
        # Below line to get elevations for subplot A.
        wshedElevs[currWshed] = collect(skipmissing(grp.Elevation))
        # Below section to get SDI across watersheds for subplot B.
        wshedSamples = unique(grp.sample)
        currWshedSDI = filter(k -> first(k) in wshedSamples, allSDIs)
        wshedSDIs[currWshed] = collect(values(currWshedSDI))
        # Below section is SDI against elevation for subplot C.
        sdiDF = DataFrame(sample = collect(keys(currWshedSDI)),
                            sdi = collect(values(currWshedSDI)))
        subSampleDB = filter(:sample => n -> in(n, wshedSamples), sampleDB)
        select!(subSampleDB, [:sample, :Elevation, :pH, :wshed])
        tmpDf = innerjoin(sdiDF, subSampleDB, on = :sample)
        plottingData = vcat(plottingData, tmpDf)
    end
    cm, pt = 96/2.54 , 4/3
    figure = Figure(fontsize = 14pt, size = (29.21cm, 35cm))
    elevData, sdiData = dictToPlotArr(wshedElevs), dictToPlotArr(wshedSDIs)

    ax1 = Axis(figure[1,1], xticks = (1:1:length(elevData.z), elevData.z), ylabel = "Elevation (m)", xticklabelrotation=45.0)
    boxplot!(ax1,elevData.x,elevData.y, strokewidth = 2)

    ax2 = Axis(figure[1,2], xticks = (1:1:length(sdiData.z), sdiData.z), ylabel = "SDI", xticklabelrotation=45.0)
    boxplot!(ax2, sdiData.x, sdiData.y, strokewidth = 2)

    ax3 = Axis(figure[2,1], ylabel = "SDI", xlabel = "Elevation (m)")
    scatter!(ax3, plottingData.Elevation, plottingData.sdi)

    ax4 = Axis(figure[2,2], xlabel = "pH")
    scatter!(ax4, plottingData.pH, plottingData.sdi)

    sampleSDIs = Dict()
    for grp in groupby(sampleDB, :SampleType)
        currSampleT = grp.SampleType[1]
        currSampleSDI = filter(k -> first(k) in unique(grp.sample), allSDIs)
        sampleSDIs[currSampleT] = collect(values(currSampleSDI))
    end
    filter!(p -> !ismissing(first(p)),  sampleSDIs)
    sampleTypeData = dictToPlotArr(sampleSDIs)

    ax5 = Axis(figure[3,1:2], xticks = (1:1:length(sampleTypeData.z), sampleTypeData.z),
               ylabel = "SDI", xticklabelrotation=45.0)
    boxplot!(ax5, sampleTypeData.x, sampleTypeData.y, strokewidth = 2)

    subSecKwargs = (fontsize = 18pt, font = :bold, halign = :right, padding = (0, 5, 5, 0))
    Label(figure[1,1, TopLeft()], "A" ; subSecKwargs...)
    Label(figure[2,1, TopLeft()], "B" ; subSecKwargs...)
    Label(figure[2,2, TopLeft()], "C" ; subSecKwargs...)
    Label(figure[3,1, TopLeft()], "D" ; subSecKwargs...)



    return figure
end
plot = main()
save(joinpath(@__DIR__, "private/", "plots", "diversity.svg"), plot)