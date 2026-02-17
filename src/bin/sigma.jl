####################################################################################################
# cli args
####################################################################################################

# Load path definitions
include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.UTIL, "args.jl"))
args = sigma_args()

###################################################################################################

include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "params.jl"))
include(joinpath(Paths.UTIL, "sigma.jl"))

###################################################################################################

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = loadTrajectoryParams(
  args["trajectory"],
  Dict(
    "metadata" => Vars.SMETA_xlsx,
    "batches" => [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
  ),
)

@vinfo args sigma_params

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

@vinfo args "Done"

####################################################################################################
