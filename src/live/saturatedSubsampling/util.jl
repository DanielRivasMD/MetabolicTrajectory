####################################################################################################

using Colors
using DataFrames
using Dates
using Plots
using Statistics
using XLSX

####################################################################################################

"""
    compare_columns(df::DataFrame, subdf::DataFrame; drop=String[])

Compare column names between the original `df` and a suffix‑stripped `subdf`.

- Strips trailing `_N` suffixes from `df` column names.
- Normalizes both sides to `Symbol` for consistent comparison.
- Optionally drops metadata columns via `drop`.

Prints the columns only in `df` (after stripping) and only in `subdf`.
"""
function compare_columns(df::DataFrame, subdf::DataFrame; drop = String[])
  # Normalize df names: strip suffixes and convert to Symbol
  df_base_syms = unique(Symbol.(replace.(String.(names(df)), r"_\d+$" => "")))
  # Drop metadata columns if requested
  df_base_syms = setdiff(df_base_syms, Symbol.(drop))

  # Normalize subdf names to Symbol
  sub_syms = Symbol.(names(subdf))

  # Differences
  only_in_df = [c for c in df_base_syms if c ∉ sub_syms]
  only_in_subdf = [c for c in sub_syms if c ∉ df_base_syms]

  return (; only_in_df, only_in_subdf)
end

"""
    compare_columns(df::DataFrame, subdfs::Dict{Int,DataFrame}; drop=String[])

Compare column names between the original `df` and each suffix‑stripped subdf in `subdfs`.
Returns a Dict keyed by animal ID with the differences.
"""
function compare_columns(df::DataFrame, subdfs::Dict{Int,DataFrame}; drop = String[])
  results = Dict{Int,NamedTuple}()
  for (id, subdf) in subdfs
    results[id] = compare_columns(df, subdf; drop = drop)
  end
  return results
end

####################################################################################################

"""
    load_experiments(p::TrajectoryParams)

Load experiments described in an XLSX file (from `p.metadata`) and
corresponding CSV batches (from `p.batches`).

Processing steps:

- **Sheet selection**: Each sheet in the XLSX is inspected; the sheet name must
  contain a date in `YYYY-MM-DD` format, which is used as the experiment key.
- **Metadata normalization**:
  - `Cage_nr` → cast to `Int`.
  - `Animal_nr` → replace `"EMPTY"` with `0`, then cast to `Int`.
  - `Sex` and `Genotype` → replace `"EMPTY"` with `""`, then cast to `String`.
- **Experiment normalization**:
  - Columns listed in `Vars.xvars_csv` are parsed as `DateTime` using the
    format `"mm/dd/yyyy HH:MM:SS"`.
  - All other columns are cast to `Float64`.
- **Matching**: The extracted date key is matched against the provided CSV
  batch filenames to locate the corresponding experiment data.
- **Output**: Returns a `Dict{String,ExperimentBundle}` keyed by the date string,
  where each `ExperimentBundle` contains the cleaned `metadata::DataFrame` and
  `experiment::DataFrame`.
"""

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

"""
    split_by_animal(bundle::ExperimentBundle; timecol::Symbol = :Date_Time)

Split the experiment DataFrame in `bundle` into per‑animal DataFrames.

- Uses suffixes in column names (e.g. `_1`, `_2`) to group signals.
- Renames columns to strip suffixes.
- Keys the output dictionary by the `Animal_nr` values from `bundle.metadata`,
  instead of the raw suffix index.
- Ignores animals with `Animal_nr == 0`.
"""
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

"""
    split_by_animal(bundles::Dict{String,ExperimentBundle}; timecol::Symbol = :Date_Time)

Split all experiment DataFrames in `bundles` into per‑animal DataFrames.

- Calls `split_by_animal(bundle)` on each entry.
- Consolidates into a single dictionary keyed by `animal_nr`.
- Ignores animals with `Animal_nr == 0`.
"""
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

"""
    collect_subsamples(subdfs::Dict{Int,DataFrame}, var::Symbol, params::TrajectoryParams)

Collect subsamples from a single variable across all animals according to parameters in `params`.

- **Reference length** is the maximum number of rows across all animals in `subdfs`.
- **`p.subsample_len`**:
  - If `0 < subsample_len < 1`, interpreted as a fraction of the reference length.
  - If `subsample_len >= 1`, interpreted as an absolute number of indices.
- **`params.subsample_var`** applies +-variance jitter to the subsample length.
- **`params.nsamples`** controls how many subsamples to draw per animal.

Each subsample is a contiguous slice of the signal, chosen at a random starting index,
with length jittered around the target subsample length.

Returns a named tuple:
- `subsamples::Vector{Vector{Float64}}` — the collected subsample vectors
- `ids::Vector{String}` — identifiers of the form `"animalID_startIndex"`.
"""
function collect_subsamples(
  subdfs::Dict{Int,DataFrame},
  var::Symbol,
  params::TrajectoryParams,
)
  # Reference length: maximum number of rows across all animals
  max_rows = maximum(nrow(df) for df in values(subdfs))

  basevar = Symbol(replace(string(var), r"_\d+$" => ""))

  all_subsamples = Vector{Vector{Float64}}()
  all_ids = String[]

  for (animal_id, subdf) in subdfs
    signal = collect(skipmissing(subdf[!, basevar]))
    n = length(signal)
    if n == 0
      @warn "Animal $animal_id has no valid data for $basevar, skipping"
      continue
    end

    # Determine subsample length
    subsample_len = if 0 < params.subsample_len < 1
      round(Int, params.subsample_len * max_rows)
    else
      round(Int, params.subsample_len)
    end

    for i = 1:params.nsamples
      # Apply variance jitter
      len_var = round(
        Int,
        subsample_len * (1 + (rand() * params.subsample_var * 2 - params.subsample_var)),
      )

      if n <= len_var
        @warn "Animal $animal_id has shorter signal ($n) than subsample length ($len_var), skipping"
        continue
      end

      start_idx = rand(1:(n-len_var))
      subsample = signal[start_idx:start_idx+len_var-1]

      push!(all_subsamples, subsample)
      push!(all_ids, string(animal_id, "_", start_idx))
    end
  end

  return (subsamples = all_subsamples, ids = all_ids)
end

"""
    collect_subsamples(subdfs::Dict{Int,DataFrame}, vars::Vector{Symbol}, params::TrajectoryParams)

Collect subsamples for multiple variables by delegating to
`collect_subsamples(subdfs, var, p)` for each variable.

Returns a dictionary keyed by normalized variable symbol (suffix `_N` stripped).
Each value is a named tuple:
- `subsamples::Vector{Vector{Float64}}`
- `ids::Vector{String}`
"""
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

"""
    cost_stats(A)

Compute basic statistics of a cost matrix (or any array):
- minimum
- maximum
- mean
- median
- standard deviation

Returns a NamedTuple.
"""
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

using Plots, Colors

"""
    plot_grouped_costmatrix(cost_matrix, groups; pad=10, group_colors=nothing)

Embed group membership as colored bars along the top and left edges of a cost matrix
and plot as a single heatmap.

Arguments:
- `cost_matrix::Matrix{Float64}` — the N×N cost matrix
- `groups::Vector{String}` — length N, group label for each trajectory
- `pad::Int` — thickness of the colored bars (default 10)
- `group_colors::Union{Nothing, Vector{<:Colorant}}` — optional custom colors for groups.
   If `nothing`, colors are generated automatically in the green/blue range.

Behavior:
- Normalizes cost values to 0–100
- Uses provided colors if given, otherwise generates them
- Renders the cost matrix and group bars in two masked layers
- Returns the plot object, group levels, and colors
"""
function plot_grouped_costmatrix(
  cost_matrix::Matrix{Float64},
  groups::Vector{String};
  pad::Int = 10,
  group_colors::Union{Nothing,Vector{<:Colorant}} = nothing,
)

  N = size(cost_matrix, 1)
  @assert length(groups) == N "groups vector must match matrix size"

  # Map groups to codes
  group_levels = unique(groups)
  group_to_code = Dict(g => i for (i, g) in enumerate(group_levels))
  codes = [group_to_code[g] for g in groups]

  # # Normalize cost matrix to 0–100
  # Cmin, Cmax = extrema(cost_matrix)
  # C_scaled = (cost_matrix .- Cmin) ./ (Cmax - Cmin) .* 100

  # Assign colors
  if group_colors === nothing
    base = RGB(0, 0.5, 0.5)
    group_colors = distinguishable_colors(length(group_levels), [base])
  else
    @assert length(group_colors) == length(group_levels) "Number of colors must match number of groups"
  end

  # --- Build masked matrices ---
  # Core only (pad = NaN)
  core_only = fill(NaN, N + pad, N + pad)
  core_only[pad+1:end, pad+1:end] .= cost_matrix

  # Pad only (core = NaN)
  pad_only = fill(NaN, N + pad, N + pad)
  for i = 1:N
    pad_only[1:pad, pad+i] .= codes[i]      # top strip
    pad_only[pad+i, 1:pad] .= codes[i]      # left strip
  end

  plt = heatmap(core_only; color = :inferno, yflip = true, legend = false)
  heatmap!(
    pad_only;
    color = cgrad(group_colors, categorical = true),
    yflip = true,
    legend = false,
  )

  return plt, group_levels, group_colors
end

using LinearAlgebra

# Classical MDS: from distance matrix D to m-dimensional embedding (default 50)
function mds_embedding(D::AbstractMatrix{<:Real}; m::Int = 50)
  @assert size(D, 1) == size(D, 2) "D must be square"
  n = size(D, 1)
  # Square distances
  D2 = D .^ 2
  # Double-centering: B = -0.5 * J * D2 * J
  J = I - (1 / n) * ones(n, n)
  B = -0.5 * J * D2 * J
  # Eigen-decomposition (symmetric)
  vals, vecs = eigen(Symmetric(B))
  # Keep positive eigenvalues (numerical safety)
  pos = findall(>(0), vals)
  if isempty(pos)
    error("MDS failed: no positive eigenvalues")
  end
  k = min(m, length(pos))
  idx = reverse(pos)[1:k]               # largest positive
  Λ = Diagonal(sqrt.(vals[idx]))        # scale by sqrt(eigenvalues)
  X = vecs[:, idx] * Λ                  # n × k embedding
  return X
end

using MultivariateStats, Clustering, Statistics

function pca_clustering_from_cost(
  D::AbstractMatrix{<:Real};
  n_components::Int = 10,
  k::Int = 5,
)
  # Embed distances into Euclidean space first
  X = mds_embedding(D; m = max(n_components, 50))

  # Center the data
  Xc = X .- mean(X, dims = 1)

  # PCA
  M = fit(PCA, Xc; maxoutdim = n_components)
  scores = MultivariateStats.transform(M, Xc)

  # Cluster
  R = kmeans(scores', k)
  return (scores = scores, assignments = R.assignments, centers = R.centers)
end

using TSne, Random, Clustering

function tsne_clustering_from_cost(
  D::AbstractMatrix{<:Real};
  k::Int = 5,
  no_dims::Int = 2,
  initial_dims::Int = 50,
  perplexity::Int = 30,
  max_iter::Int = 1000,
  seed::Int = 42,
)

  # Embed distances into Euclidean space first
  X = mds_embedding(D; m = initial_dims)

  # Run t-SNE with perplexity as positional argument
  Random.seed!(seed)
  Y = tsne(X, no_dims, initial_dims, perplexity, max_iter)

  # Cluster the embedding
  R = kmeans(Y', k)
  return (embedding = Y, assignments = R.assignments)
end

###################################################################################################
