using CSV, DataFrames, CairoMakie, StatsBase


binPath = joinpath(@__DIR__,"data", "outputDB.csv")
binDB = coalesce.(CSV.read(binPath, DataFrame), "Missing")
