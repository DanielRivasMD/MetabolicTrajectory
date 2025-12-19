# Load path definitions
include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.SRC, "live/saturatedSubsampling/config.jl"))
include(joinpath(Paths.SRC, "live/saturatedSubsampling/util.jl"))

# Load configuration structs
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))

using Clustering
using Dates
using DataFrames
using Distances
using DynamicAxisWarping
using LinearAlgebra
using Plots
using Statistics
using Random
using UMAP

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = TrajectoryParams(
  metadata = Vars.SMETA_xlsx,
  batches = [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
)

# Load and split
bundles = load_experiments(sigma_params)

# collect metadata
# Collect all metadata into one DataFrame
meta = vcat([b.metadata for b in values(bundles)]...)

# Keep only the columns of interest
meta = select(meta, [:Animal_nr, :Sex, :Genotype])

# Drop rows where Animal_nr == 0 (if those are placeholders)
filter!(row -> row.Animal_nr != 0, meta)
rename!(meta, :Animal_nr => :Animal)
meta.Group = string.(meta.Sex, "_", meta.Genotype)
animal_to_group = Dict(row.Animal => row.Group for row in eachrow(meta))

subdfs = split_by_animal(bundles)

using DataFrames, Dates

function summarize_first4h(subdfs::Dict{Int,DataFrame})
  results = DataFrame(ID = Int[], RT_PedMeters = Float64[], RT_Water = Float64[])
  for (id, df) in subdfs
    # find the time window
    start_time = first(df.Date_Time)
    cutoff = start_time + Hour(4)

    # slice the first 4 hours
    df4h = filter(row -> row.Date_Time <= cutoff, df)

    # cumulative sums
    ped_sum = sum(df4h.RT_PedMeters)
    water_sum = sum(df4h.RT_Water)

    push!(results, (id, ped_sum, water_sum))
  end
  return results
end

summary_df = summarize_first4h(subdfs)

summary = leftjoin(summary_df, meta[:, [:Animal, :Group]], on = :ID => :Animal)

