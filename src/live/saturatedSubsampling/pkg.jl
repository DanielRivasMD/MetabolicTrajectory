# Load path definitions
include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

# Load configuration structs
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))

using Clustering
using Dates
using DataFrames
using Distances
using DynamicAxisWarping
using Plots
using Statistics
using Random


# Load data
df = load_timeseries(Vars.SIG1R_HT_file)

# Group by subject
gdf = groupby(df, :Animal)

# Variables
vars = setdiff(names(df), Vars.exclude_vars)
