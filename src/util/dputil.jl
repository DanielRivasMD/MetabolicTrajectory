####################################################################################################

module DPCore

####################################################################################################

using Dates
using DataFrames
using XLSX
using TOML
using ..USCore: readdf, writedf

####################################################################################################

const TIME_COLS = ["Date_Time"]

####################################################################################################

struct IngestParams
  metadata::String
  batches::Vector{String}
  out_dir::String
end

####################################################################################################

struct ExperimentBundle
  metadata::DataFrame
  experiment::DataFrame
end

####################################################################################################

default_out_dir() = "sigma"

####################################################################################################

function load_ingest_params(config_path::String)
  raw = TOML.parsefile(config_path)
  section = get(raw, "ingest", raw)
  metadata = section["metadata"]
  batches = section["batches"]
  out_dir = get(section, "out_dir", default_out_dir())
  return IngestParams(metadata, batches, out_dir)
end

####################################################################################################

"""
    run_ingestion(params::IngestParams) -> Dict

Performs the full ingestion:
  1. Load experiments from metadata XLSX and CSV batches.
  2. Build a unified metadata DataFrame (Animal, Sex, Genotype).
  3. Split time‑series data by animal.
  4. Write meta.csv and per‑animal CSVs to `params.out_dir`.

Returns a dictionary with `out_dir` and `animal_count`.
"""
function run_ingestion(params::IngestParams)
  mkpath(params.out_dir)

  bundles = load_experiments(params)
  meta = build_metadata(bundles)

  # Save metadata
  writedf(joinpath(params.out_dir, "meta.csv"), meta; sep = ',')

  # Split and write per‑animal CSVs
  dfs = split_by_animal(bundles)
  writedf_dict(params.out_dir, dfs; sep = ',')

  return Dict("out_dir" => params.out_dir, "animal_count" => length(dfs))
end

####################################################################################################

function build_metadata(bundles::Dict{String,ExperimentBundle})
  meta = vcat([b.metadata for b in values(bundles)]...)
  meta = select(meta, [:Animal_nr, :Sex, :Genotype])
  filter!(row -> row.Animal_nr != 0, meta)
  rename!(meta, :Animal_nr => :Animal)
  meta.Group = string.(meta.Sex, "_", meta.Genotype)
  return meta
end

####################################################################################################

function load_experiments(params::IngestParams)
  return load_experiments(params.metadata, params.batches)
end

####################################################################################################

function load_experiments(metadata_path::String, batches::Vector{String})
  bundles = Dict{String,ExperimentBundle}()

  XLSX.openxlsx(metadata_path) do xf
    for sheetname in XLSX.sheetnames(xf)
      m = match(r"\d{4}-\d{2}-\d{2}", sheetname)
      isnothing(m) && continue
      datekey = m.match

      sheet = xf[sheetname]
      metadata = XLSX.gettable(sheet) |> DataFrame

      # Normalise metadata columns
      rename!(metadata, Symbol.(names(metadata)))
      metadata[!, :Cage_nr] = Int.(metadata[!, :Cage_nr])
      metadata[!, :Animal_nr] =
        map(x -> x == "EMPTY" ? 0 : parse(Int, string(x)), metadata[!, :Animal_nr])
      metadata[!, :Sex] = map(x -> x == "EMPTY" ? "" : string(x), metadata[!, :Sex])
      metadata[!, :Genotype] =
        map(x -> x == "EMPTY" ? "" : string(x), metadata[!, :Genotype])

      # Find matching CSV
      csvmatch = filter(f -> occursin(datekey, f), batches)
      if isempty(csvmatch)
        @warn "No CSV file found for date $datekey"
        continue
      end

      df = readdf(first(csvmatch); sep = ',')

      vars = setdiff(names(df), TIME_COLS)

      # Parse time columns
      for c in TIME_COLS
        df[!, c] = DateTime.(df[!, c], dateformat"mm/dd/yyyy HH:MM:SS")
      end

      # Cast all other columns to Float64
      for c in vars
        df[!, c] = Float64.(df[!, c])
      end

      bundles[datekey] = ExperimentBundle(metadata, df)
    end
  end

  return bundles
end

####################################################################################################

function split_by_animal(bundle::ExperimentBundle; timecol::Symbol = Symbol(TIME_COLS[1]))
  df = bundle.experiment
  meta = bundle.metadata

  # Map column suffixes to animal IDs
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
    if animal_nr == 0
      continue
    end

    sdf = df[:, vcat([timecol], sort(cols))]
    # Strip suffix from column names
    renames = Dict(c => Symbol(replace(String(c), r"_\d+$" => "")) for c in cols)
    rename!(sdf, renames)
    bases = sort(Symbol.(replace.(String.(cols), r"_\d+$" => "")))
    subdfs[animal_nr] = sdf[:, vcat([timecol], bases)]
  end

  return subdfs
end

####################################################################################################

function split_by_animal(
  bundles::Dict{String,ExperimentBundle};
  timecol = Symbol(TIME_COLS[1]),
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

####################################################################################################

function writedf_dict(dir::String, dict::Dict{Int,DataFrame}; sep = ',')
  for (k, df) in dict
    path = joinpath(dir, "exp_$(k).csv")
    writedf(path, df; sep = sep)
  end
end

####################################################################################################

function load_timeseries(path::AbstractString)
  df = readxlsx(path, sheetname = "TimeSeries")
  df.DateTime = DateTime.(df.DateTime, dateformat"yyyy/mm/dd HH:MM:SS")
  for c in names(df)
    col = df[!, c]
    if c == :Animal
      df[!, c] = string.(col)
    elseif eltype(col) <: AbstractString
      df[!, c] = map(x -> x == "." ? missing : tryparse(Float64, x), col)
    elseif eltype(col) == Any
      df[!, c] = map(x -> x == "." ? missing : tryparse(Float64, string(x)), col)
    end
  end
  return df
end

####################################################################################################

end

####################################################################################################
