# # interactive
# include("src/live/saturatedSubsampling/config.jl")
# include("src/live/saturatedSubsampling/util.jl")

# Load path definitions
include(joinpath((pwd() == @__DIR__) ? "src" : "../..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

# Load configuration structs
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))

@info "paths loaded"

using Clustering
using Dates
using DataFrames
using Distances
using DynamicAxisWarping
using Plots
using Random

@info "dependencies loaded"

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = TrajectoryParams(
  metadata = Vars.SMETA_xlsx,
  batches = [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
)

# Load and split
bundles = load_experiments(sigma_params)

# Collect all metadata into one DataFrame
meta = vcat([b.metadata for b in values(bundles)]...)

# Keep only the columns of interest
meta = select(meta, [:Animal_nr, :Sex, :Genotype])

# Drop rows where Animal_nr == 0 (if those are placeholders)
filter!(row -> row.Animal_nr != 0, meta)
rename!(meta, :Animal_nr => :Animal)

subdfs = split_by_animal(bundles)

# Choose variables
vars = names(subdfs[1])[2:end] .|> Symbol

for var in vars

  cost_path = joinpath(Paths.TMP, string(var, "_cost_matrix.csv"))
  ids_path = joinpath(Paths.TMP, string(var, "_order_ids.csv"))

  @info "Loading cached DTW results for $var"

  cost_matrix = readdlm(cost_path, ',', Float64)
  all_ordered_ids = dataframe_to_ids(readdf(ids_path; sep = ','))

  plt = plot_grouped_costmatrix(cost_matrix, all_ordered_ids, meta; pad = 200, title = String(var))
  display(plt)

  # Hierarchical clustering
  tree = hclust(cost_matrix; linkage = :ward)

  plt = plot_grouped_costmatrix(
  cost_matrix[tree.order, tree.order],
    all_ordered_ids[tree.order],
    meta;
    pad = 200,
    title = String(var),
  )
  display(plt)

end
