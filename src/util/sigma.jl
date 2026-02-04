####################################################################################################

"""
    ExperimentBundle(metadata::DataFrame, experiment::DataFrame)

Container for a single experiment, holding:
- `metadata`  : DataFrame describing cages, animals, sex, and group (from XLSX sheet)
- `experiment`: DataFrame with the corresponding experimental time‑series data (from CSV)
"""
struct ExperimentBundle
  metadata::DataFrame
  experiment::DataFrame
end

####################################################################################################

function load_experiments(params::TrajectoryParams)
  bundles = Dict{String,ExperimentBundle}()

  XLSX.openxlsx(params.metadata) do xf
    for sheetname in XLSX.sheetnames(xf)
      # Extract date in YYYY-MM-DD format from sheet name
      m = match(r"\d{4}-\d{2}-\d{2}", sheetname)
      isnothing(m) && continue
      datekey = m.match

      # Load and normalize metadata
      sheet = xf[sheetname]
      metadata = XLSX.gettable(sheet) |> DataFrame

      # Ensure expected columns exist
      rename!(metadata, Symbol.(names(metadata)))  # force Symbol names

      # Cast Cage_nr to Int
      metadata[!, :Cage_nr] = Int.(metadata[!, :Cage_nr])

      # Cast Animal_nr: replace "EMPTY" with 0, then Int
      metadata[!, :Animal_nr] =
        map(x -> x == "EMPTY" ? 0 : parse(Int, string(x)), metadata[!, :Animal_nr])

      # Cast Sex: replace "EMPTY" with "" and String
      metadata[!, :Sex] = map(x -> x == "EMPTY" ? "" : string(x), metadata[!, :Sex])

      # Cast Genotype: replace "EMPTY" with "" and String
      metadata[!, :Genotype] =
        map(x -> x == "EMPTY" ? "" : string(x), metadata[!, :Genotype])

      # Find matching CSV file
      csvmatch = filter(f -> occursin(datekey, f), params.batches)
      if isempty(csvmatch)
        @warn "No CSV file found for date $datekey"
        continue
      end

      # Load and normalize experiment
      df = readdf(first(csvmatch); sep = ',')

      vars = setdiff(names(df), Vars.xvars_csv)

      # Parse time columns
      for c in Vars.xvars_csv
        df[!, c] = DateTime.(df[!, c], dateformat"mm/dd/yyyy HH:MM:SS")
      end

      # Cast all other columns to Float64
      for c in vars
        df[!, c] = Float64.(df[!, c])
      end

      # Store bundle
      bundles[datekey] = ExperimentBundle(metadata, df)
    end
  end

  return bundles
end

####################################################################################################

function split_by_animal(bundle::ExperimentBundle; timecol::Symbol = :Date_Time)
  df = bundle.experiment
  meta = bundle.metadata

  # Map columns to suffix IDs
  suffix_ids = Dict{Symbol,Int}()
  for col in names(df)
    m = match(r"_(\d+)$", String(col))
    if m !== nothing
      suffix_ids[Symbol(col)] = parse(Int, m.captures[1])
    end
  end

  grouped = Dict{Int,Vector{Symbol}}()
  for (col, id) in suffix_ids
    push!(get!(grouped, id, Symbol[]), col)
  end

  # Build suffix → Animal_nr mapping
  raw_ids = meta[!, :Animal_nr]
  idmap = Dict(i => raw_ids[i] for i = 1:nrow(meta))

  subdfs = Dict{Int,DataFrame}()
  for (suffix_id, cols) in grouped
    animal_nr = idmap[suffix_id]
    # Skip empty/ignored animals
    if animal_nr == 0
      continue
    end

    sdf = df[:, vcat([timecol], sort(cols))]
    renames = Dict(c => Symbol(replace(String(c), r"_\d+$" => "")) for c in cols)
    rename!(sdf, renames)
    bases = sort(Symbol.(replace.(String.(cols), r"_\d+$" => "")))
    subdfs[animal_nr] = sdf[:, vcat([timecol], bases)]
  end

  return subdfs
end

function split_by_animal(
  bundles::Dict{String,ExperimentBundle};
  timecol::Symbol = :Date_Time,
)
  consolidated = Dict{Int,DataFrame}()

  for bundle in values(bundles)
    subdfs = split_by_animal(bundle; timecol = timecol)
    for (animal_nr, sdf) in subdfs
      if haskey(consolidated, animal_nr)
        @warn "Duplicate Animal_nr $animal_nr across bundles; overwriting"
      end
      consolidated[animal_nr] = sdf
    end
  end

  return consolidated
end

###################################################################################################

function collect_subsamples(signal::Vector{Float64}, params::TrajectoryParams)
  n = length(signal)
  if n == 0
    return (subsamples = Vector{Vector{Float64}}(), ids = String[])
  end

  # Determine base subsample length
  len = if 0 < params.len < 1
    round(Int, params.len * n)
  else
    round(Int, params.len)
  end

  subsamples = Vector{Vector{Float64}}()
  ids = String[]

  for i = 1:params.nsamples
    # Apply variance jitter
    len_var = round(
      Int,
      len * (1 + (rand() * params.var * 2 - params.var)),
    )

    if n <= len_var
      @warn "Signal shorter ($n) than subsample length ($len_var), skipping"
      continue
    end

    start_idx = rand(1:(n-len_var))
    subsample = signal[start_idx:start_idx+len_var-1]

    push!(subsamples, subsample)
    push!(ids, string(start_idx))
  end

  return (subsamples = subsamples, ids = ids)
end

function collect_subsamples(
  subdfs::Dict{Int,DataFrame},
  var::Symbol,
  params::TrajectoryParams,
)
  max_rows = maximum(nrow(df) for df in values(subdfs))
  basevar = Symbol(replace(string(var), r"_\d+$" => ""))

  all_subsamples = Vector{Vector{Float64}}()
  all_ids = String[]

  for (animal_id, subdf) in subdfs
    signal = collect(skipmissing(subdf[!, basevar]))
    if isempty(signal)
      @warn "Animal $animal_id has no valid data for $basevar, skipping"
      continue
    end

    # Use the vector-based function
    result = collect_subsamples(signal, params)
    append!(all_subsamples, result.subsamples)
    # Prefix IDs with animal_id
    append!(all_ids, [string(animal_id, "_", id) for id in result.ids])
  end

  return (subsamples = all_subsamples, ids = all_ids)
end

function collect_subsamples(
  subdfs::Dict{Int,DataFrame},
  vars::Vector{Symbol},
  params::TrajectoryParams,
)
  Dict(var => collect_subsamples(subdfs, var, params) for var in vars)
end

###################################################################################################

function split_id(id::String)
  parts = split(id, "_")
  return (parse(Int, parts[1]), parse(Int, parts[2]))
end

###################################################################################################

function cost_stats(A)
  return (
    min = minimum(A),
    max = maximum(A),
    mean = mean(A),
    median = median(A),
    std = std(A),
  )
end

###################################################################################################

function plot_grouped_costmatrix(
  cost_matrix::Matrix{Float64},
  groups::Vector{String};
  pad::Int = 10,
  group_colors::Union{Nothing,Vector{<:Colorant}} = nothing,
)
  N = size(cost_matrix, 1)
  @assert size(cost_matrix, 2) == N "cost_matrix must be square"
  @assert length(groups) == N "groups vector must match matrix size"

  # Map groups to integer codes
  group_levels = unique(groups)
  group_to_code = Dict(g => i for (i, g) in enumerate(group_levels))
  codes = [group_to_code[g] for g in groups]

  # Colors per group
  if group_colors === nothing
    base = RGB(0.0, 0.5, 0.5)
    group_colors = distinguishable_colors(length(group_levels), [base])
  else
    @assert length(group_colors) == length(group_levels)
  end

  # Build masked matrices
  core_only = fill(Float64(NaN), N + pad, N + pad)
  core_only[pad+1:end, 1:N] .= cost_matrix   # core in bottom-left block

  pad_colors = fill(RGBA(0, 0, 0, 0), N + pad, N + pad)
  for i = 1:N
    pad_colors[1:pad, i] .= group_colors[codes[i]]        # top strip
    pad_colors[pad+i, N+1:N+pad] .= group_colors[codes[i]] # right strip
  end

  # Plot core heatmap with axes disabled
  plt = heatmap(
    core_only;
    color = :inferno,
    yflip = true,
    legend = false,
    nan_color = RGBA(0, 0, 0, 0),
    axis = false,   # remove x/y axes annotations
  )

  # Overlay pad colors directly
  plot!(plt, pad_colors; seriestype = :heatmap, yflip = true, legend = false, axis = false)

  return plt, group_levels, group_colors
end

###################################################################################################
