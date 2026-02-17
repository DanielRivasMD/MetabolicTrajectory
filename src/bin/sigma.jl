####################################################################################################
#######################################       MODULE       #########################################
####################################################################################################

module Sigma

export run, TrajectoryParams, loadTrajectoryParams

####################################################################################################

include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.UTIL, "args.jl"))
args = sigma_args()

###################################################################################################

include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "params.jl"))
include(joinpath(Paths.UTIL, "wrangle.jl"))

###################################################################################################

# TODO: divide params for sigma & params for saturated. read from a single file, but with different sections
function run(_::Dict)

  # Sigma experiment: one metadata XLSX, three CSV batches
  sigma_params = loadTrajectoryParams(
    args["config"],
    Dict(
      "metadata" => Vars.SMETA_xlsx,
      "batches" => [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
    ),
  )

  isdir(sigma_params.meta_path) || mkpath(sigma_params.meta_path)

  ###################################################################################################

  # Load and split
  bundles = load_experiments(sigma_params)

  # Collect all metadata into one DataFrame
  meta = vcat([b.metadata for b in values(bundles)]...)

  # Keep only the columns of interest
  meta = select(meta, [:Animal_nr, :Sex, :Genotype])

  # Drop rows where Animal_nr == 0 (if those are placeholders)
  filter!(row -> row.Animal_nr != 0, meta)
  rename!(meta, :Animal_nr => :Animal)

  dfs = split_by_animal(bundles)

  writedf(joinpath(sigma_params.meta_path, "meta.csv"), meta)
  writedf_dict(sigma_params.meta_path, dfs)

end

end

####################################################################################################
#######################################       MODULE       #########################################
####################################################################################################

if abspath(PROGRAM_FILE) == @__FILE__
  Sigma.run(Dict())
end

####################################################################################################
