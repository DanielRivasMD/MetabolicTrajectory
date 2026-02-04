####################################################################################################

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = TrajectoryParams(
  metadata = Vars.SMETA_xlsx,
  batches = [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
)

@info sigma_params

sigma_params = loadTparams(args["trajectory"], sigma_params)

@info sigma_params

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

# Collect subsamples for all variables at once
subsample_results = collect_subsamples(subdfs, vars, sigma_params)

####################################################################################################
