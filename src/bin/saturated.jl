####################################################################################################
# cli args
####################################################################################################

# Load path definitions
include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.UTIL, "args.jl"))
args = saturated_args()

###################################################################################################

include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "params.jl"))
include(joinpath(Paths.UTIL, "wrangle.jl"))

####################################################################################################

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = loadTrajectoryParams(
  args["config"],
  Dict(
    "metadata" => Vars.SMETA_xlsx,
    "batches" => [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
  ),
)

@vinfo args sigma_params

####################################################################################################

dfs = readdf_dict(sigma_params.out)

# Collect subsamples for all variables at once
subsamples = collect_subsamples(dfs, vars, sigma_params)

tag = matrix_tag(sigma_params)
outdir = joinpath(Paths.TMP, tag)

# Create directory if missing
isdir(outdir) || mkpath(outdir)

for var in vars

  @vinfo args "Processing variable: $var"

  order = sortperm(subsamples[var].ids; by = id -> (id.subject, id.ixs[1]))
  ordered_subsamples = subsamples[var].subsamples[order]
  ordered_ids = subsamples[var].ids[order]

  N = length(ordered_subsamples)
  if N == 0
    @warn "No subsamples collected for $var, skipping"
  end

  # Build paths using tag directory
  cost_path = joinpath(outdir, string(var, "_cost_matrix.csv"))
  ids_path = joinpath(outdir, string(var, "_order_ids.csv"))

  @vinfo args "Computing DTW cost matrix for $var"

  # Compute DTW cost matrix
  cost_matrix = zeros(Float64, N, N)

  for i = 1:N
    for j = i:N
      cost, _, _ = dtw(ordered_subsamples[i], ordered_subsamples[j], SqEuclidean())
      norm_cost =
        cost / mean([length(ordered_subsamples[i]), length(ordered_subsamples[j])])
      cost_matrix[i, j] = norm_cost
      cost_matrix[j, i] = norm_cost
    end
  end

  # Save results
  writedlm(cost_path, cost_matrix, ',')
  writedf(ids_path, ids_to_dataframe(ordered_ids))

end

@vinfo args "Files written: $(sigma_params.out)"

####################################################################################################
