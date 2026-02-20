####################################################################################################

using Dates
using DataFrames
using Printf

####################################################################################################

struct SubSampleID
  subject::Int
  ixs::Tuple{Int,Int}
  time::Tuple{DateTime,DateTime}
end

struct SubSampleContainer
  subsamples::Vector{Vector{Float64}}
  ids::Vector{SubSampleID}
end

####################################################################################################

function merge_subsamplecontainers(
  containers::Vector{SubSampleContainer},
  subjects::Vector{Int},
)
  all_subs = Vector{Vector{Float64}}()
  all_ids = Vector{SubSampleID}()

  for (i, container) in enumerate(containers)
    subject = subjects[i]

    for (sub, id) in zip(container.subsamples, container.ids)
      # overwrite subject ID
      new_id = SubSampleID(subject, id.ixs, id.time)
      push!(all_subs, sub)
      push!(all_ids, new_id)
    end
  end

  return SubSampleContainer(all_subs, all_ids)
end

####################################################################################################

function ids_to_dataframe(ids::Vector{SubSampleID})
  return DataFrame(
    subject = getfield.(ids, :subject),
    start_idx = first.(getfield.(ids, :ixs)),
    end_idx = last.(getfield.(ids, :ixs)),
    start_time = first.(getfield.(ids, :time)),
    end_time = last.(getfield.(ids, :time)),
  )
end

function dataframe_to_ids(df::DataFrame)
  ids = SubSampleID[]
  for row in eachrow(df)
    push!(
      ids,
      SubSampleID(
        Int(row.subject),
        (Int(row.start_idx), Int(row.end_idx)),
        (DateTime(row.start_time), DateTime(row.end_time)),
      ),
    )
  end
  return ids
end

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

function collect_subsamples(
  signal::Vector{Float64},
  times::Vector{DateTime},
  params::TrajectoryParams,
)
  if params.limits != (0, 0)
    signal = signal[params.limits[1]:params.limits[2]]
  end
  n = length(signal)
  if n == 0
    return SubSampleContainer(Vector{Float64}[], SubSampleID[])
  end

  # Determine base subsample length
  len = if 0 < params.len < 1
    round(Int, params.len * n)
  else
    round(Int, params.len)
  end

  subsamples = Vector{Vector{Float64}}()
  ids = Vector{SubSampleID}()

  for _ = 1:params.nsamples
    # jitter
    len_var = round(Int, len * (1 + (rand() * params.var * 2 - params.var)))

    if n <= len_var
      @warn "Signal shorter ($n) than subsample length ($len_var), skipping"
      continue
    end

    start_idx = rand(1:(n-len_var))
    end_idx = start_idx + len_var - 1

    push!(subsamples, signal[start_idx:end_idx])

    push!(ids, SubSampleID(-1, (start_idx, end_idx), (times[start_idx], times[end_idx])))
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
  ids::Vector{SubSampleID},
  meta::DataFrame;
  pad::Int = 20,
  title::String = "",
)
  N = size(cost_matrix, 1)
  @assert size(cost_matrix, 2) == N
  @assert length(ids) == N

  # metadata lookup
  sex_lookup = Dict(row.Animal => row.Sex for row in eachrow(meta))
  genotype_lookup = Dict(row.Animal => row.Genotype for row in eachrow(meta))

  # Gradient from normalized midpoint
  start_idxs = first.(getfield.(ids, :ixs))
  end_idxs = last.(getfield.(ids, :ixs))
  midpoints = (start_idxs .+ end_idxs) ./ 2
  max_end = maximum(end_idxs)
  gradient = round.(Int, 100 .* midpoints ./ max_end)

  # Color logic
  sex_color(sex) =
    sex == "F" ? RGB(1, 0, 0) : sex == "M" ? RGB(0, 0, 1) : RGB(0.5, 0.5, 0.5)

  geno_color(g) =
    g == "S1RKO" ? RGB(0, 0, 0) : g == "WT" ? RGB(1, 1, 1) : RGB(0.5, 0.5, 0.5)

  gradient_color(v) = begin
    x = clamp(v % 100, 1, 100)
    t = x / 100
    RGB(0.2 * (1 - t), 0.8 * t, 0.2 * (1 - t))
  end

  # Time cycle logic
  function time_cycle_color(tstart::DateTime, tend::DateTime)
    is_day(h) = 6 <= h < 18
    midpoint = tstart + (tend - tstart) ÷ 2
    mid_is_day = is_day(hour(midpoint))
    mid_is_day ? RGB(1.0, 1.0, 0.0) : RGB(0.6, 0.6, 0.0)
  end

  # padded matrices
  core_only = fill(Float64(NaN), N + pad, N + pad)
  core_only[pad+1:end, 1:N] .= cost_matrix

  pad_colors = fill(RGBA(0, 0, 0, 0), N + pad, N + pad)

  # Dynamic section indexing
  # pad is divided into 4 equal blocks
  block = pad ÷ 4
  @assert block ≥ 1 "pad must be at least 4"

  sec1 = 1:block
  sec2 = block+1:2block
  sec3 = 2block+1:3block
  sec4 = 3block+1:4block   # up to pad

  for i = 1:N
    subject = ids[i].subject
    sex = sex_lookup[subject]
    geno = genotype_lookup[subject]
    v = gradient[i]

    tstart, tend = ids[i].time
    tcolor = time_cycle_color(tstart, tend)

    # Section 1: Sex
    pad_colors[sec1, i] .= sex_color(sex)
    pad_colors[pad+i, N.+sec1] .= sex_color(sex)

    # Section 2: Genotype
    pad_colors[sec2, i] .= geno_color(geno)
    pad_colors[pad+i, N.+sec2] .= geno_color(geno)

    # Section 3: Gradient
    pad_colors[sec3, i] .= gradient_color(v)
    pad_colors[pad+i, N.+sec3] .= gradient_color(v)

    # Section 4: Time cycle
    pad_colors[sec4, i] .= tcolor
    pad_colors[pad+i, N.+sec4] .= tcolor
  end

  plt = heatmap(
    core_only;
    title = title,
    color = :inferno,
    yflip = true,
    size = (800, 700),
    legend = false,
    nan_color = RGBA(0, 0, 0, 0),
    axis = false,
  )

  plot!(plt, pad_colors; seriestype = :heatmap, yflip = true, legend = false, axis = false)

  return plt
end

###################################################################################################

function matrix_tag(p::TrajectoryParams)
  return @sprintf(
    "n%d_len%.3f_var%.3f_limits%d_%d",
    p.nsamples,
    p.len,
    p.var,
    p.limits[1],
    p.limits[2]
  )
end

###################################################################################################
