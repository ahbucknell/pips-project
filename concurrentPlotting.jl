using CSV, CairoMakie, DataFrames

mxPath = joinpath(@__DIR__, "private", "data", "concurrencyMatrix.csv")
mx = Matrix{Int}(CSV.read(mxPath, DataFrame))

heatmap(mx)

m