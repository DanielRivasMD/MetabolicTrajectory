####################################################################################################

# TODO: to be tested
# Load path definitions
include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

####################################################################################################

include(joinpath(Paths.CONFIG, "dependencies.jl"))
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.CONFIG, "tparams.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "sigma.jl"))

####################################################################################################

include(joinpath(Paths.PIPE, "sigma.jl"))

####################################################################################################

for var in vars

  cost_path = joinpath(Paths.TMP, string(var, "_cost_matrix.csv"))
  ids_path = joinpath(Paths.TMP, string(var, "_order_ids.csv"))

  @info "Loading cached DTW results for $var"

  cost_matrix = readdlm(cost_path, ',', Float64)
  all_ordered_ids = dataframe_to_ids(readdf(ids_path; sep = ','))

  plt = plot_grouped_costmatrix(
    cost_matrix,
    all_ordered_ids,
    meta;
    pad = 200,
    title = String(var),
  )
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

####################################################################################################
