####################################################################################################

using Colors
using DataFrames
using Dates
using Plots
using Statistics
using XLSX

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

function collect_subsamples(
  signal::Vector{Float64},
  times::Vector{DateTime},
  params::TrajectoryParams,
)
  n = length(signal)
  if n == 0
    return SubSampleContainer(Vector{Float64}[], SubSampleID[])
  end

  # Determine base subsample length
  subsample_len = if 0 < params.subsample_len < 1
    round(Int, params.subsample_len * n)
  else
    round(Int, params.subsample_len)
  end

  subsamples = Vector{Vector{Float64}}()
  ids = Vector{SubSampleID}()

  for _ = 1:params.nsamples
    # jitter
    len_var = round(
      Int,
      subsample_len * (1 + (rand() * params.subsample_var * 2 - params.subsample_var)),
    )

    if n <= len_var
      @warn "Signal shorter ($n) than subsample length ($len_var), skipping"
      continue
    end

    start_idx = rand(1:(n-len_var))
    end_idx = start_idx + len_var - 1

    push!(subsamples, signal[start_idx:end_idx])

    push!(ids, SubSampleID(
      -1,
      (start_idx, end_idx),
      (times[start_idx], times[end_idx]),
    ))
  end

  return SubSampleContainer(subsamples, ids)
end

function collect_subsamples(
  subdfs::Dict{Int,DataFrame},
  var::Symbol,
  params::TrajectoryParams,
)
  basevar = Symbol(replace(string(var), r"_\d+$" => ""))

  all_subsamples = Vector{Vector{Float64}}()
  all_ids = Vector{SubSampleID}()

  for (animal_id, subdf) in subdfs
    signal = collect(skipmissing(subdf[!, basevar]))
    times = collect(skipmissing(subdf[!, :Date_Time]))

    if isempty(signal)
      @warn "Animal $animal_id has no valid data for $basevar, skipping"
      continue
    end

    result = collect_subsamples(signal, times, params)

    append!(all_subsamples, result.subsamples)

    # update subject ID
    for id in result.ids
      push!(all_ids, SubSampleID(Int32(animal_id), id.ixs, id.time))
    end
  end

  return SubSampleContainer(all_subsamples, all_ids)
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

function group_pairs(v::Vector{String})
  d = Dict{Int,Vector{Int}}()

  for s in v
    parts = split(s, "_")
    @assert length(parts) == 2 "String '$s' must contain exactly one '_'"

    key = parse(Int, parts[1])
    val = parse(Int, parts[2])

    if !haskey(d, key)
      d[key] = Int[]
    end

    push!(d[key], val)
  end

  return d
end

using Plots

function histograms_from_dict(d::Dict{Int,Vector{Int}})
  keys_sorted = sort(collect(keys(d)))

  # one subplot per key
  n = length(keys_sorted)
  plt = plot(layout = (n, 1), size = (600, 200 * n))

  for (i, k) in enumerate(keys_sorted)
    vals = d[k]

    histogram!(
      plt[i],
      vals;
      bins = :auto,
      title = "Key = $k",
      xlabel = "Values",
      ylabel = "Count",
      legend = false,
      color = :steelblue,
    )
  end

  return plt
end

###################################################################################################

using Plots
using Colors

function plot_grouped_costmatrix(
  cost_matrix::Matrix{Float64},
  groups::Vector{String},
  gradient::Vector{Int64};
  pad::Int = 15,
)
  N = size(cost_matrix, 1)
  @assert size(cost_matrix, 2) == N "cost_matrix must be square"
  @assert length(groups) == N "groups vector must match matrix size"
  @assert length(gradient) == N "gradient vector must match matrix size"

  # --- Hardcoded color logic ---

  # Section 1: Sex (F/M)
  sex_color(g) =
    occursin("F", g) ? RGB(1, 0, 0) : occursin("M", g) ? RGB(0, 0, 1) : RGB(0.5, 0.5, 0.5)

  # Section 2: Genotype (S1RKO / WT)
  geno_color(g) =
    occursin("S1RKO", g) ? RGB(0, 0, 0) :
    occursin("WT", g) ? RGB(1, 1, 1) : RGB(0.5, 0.5, 0.5)

  # Section 3: Gradient (1–100)
  gradient_color(v) = begin
    x = v % 100                     # last two digits
    x = clamp(x, 1, 100)            # ensure 1–100
    t = x / 100                     # normalize 0–1
    RGB(0.2 * (1-t), 0.8 * t, 0.2 * (1-t))  # light green -> dark green
  end

  # --- Build padded matrices ---
  core_only = fill(Float64(NaN), N + pad, N + pad)
  core_only[pad+1:end, 1:N] .= cost_matrix

  pad_colors = fill(RGBA(0, 0, 0, 0), N + pad, N + pad)

  # Three 5‑pixel sections
  sec1 = 1:5
  sec2 = 6:10
  sec3 = 11:15

  for i = 1:N
    g = groups[i]
    v = gradient[i]

    # Section 1: Sex
    pad_colors[sec1, i] .= sex_color(g)
    pad_colors[pad+i, N.+sec1] .= sex_color(g)

    # Section 2: Genotype
    pad_colors[sec2, i] .= geno_color(g)
    pad_colors[pad+i, N.+sec2] .= geno_color(g)

    # Section 3: Gradient
    pad_colors[sec3, i] .= gradient_color(v)
    pad_colors[pad+i, N.+sec3] .= gradient_color(v)
  end

  # Plot core heatmap
  plt = heatmap(
    core_only;
    color = :inferno,
    yflip = true,
    size = (800, 700),
    legend = false,
    nan_color = RGBA(0, 0, 0, 0),
    axis = false,
  )

  # Overlay padding
  plot!(plt, pad_colors; seriestype = :heatmap, yflip = true, legend = false, axis = false)

  return plt
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
