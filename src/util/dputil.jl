####################################################################################################

module DPCore

####################################################################################################

using Dates
using DataFrames
using TOML
using XLSX
using ..USCore: readdf, readxlsx, writedf

####################################################################################################

const TIME_COLS = ["Date_Time"]
const DEFAULT_OUTMETA = "meta.csv"
const ANIMAL_NR = Symbol("Animal_nr")

####################################################################################################

struct ExperimentBundle
  metadata::DataFrame
  experiment::DataFrame
end

####################################################################################################

struct DPParams
  expmeta::String
  batches::Vector{String}
  outdir::String
  outmeta::String
end

####################################################################################################

function load_params(config_path::String)
  raw = TOML.parsefile(config_path)
  # TODO: check that fields are present, if not error out with docs
  section = get(raw, "data_processing", raw)
  expmeta = section["expmeta"]
  batches = section["batches"]
  outdir = section["outdir"]
  outmeta = get(section, "outmeta", DEFAULT_OUTMETA)
  return DPParams(expmeta, batches, outdir, outmeta)
end

####################################################################################################

function build_meta(bundles::Dict{String,ExperimentBundle})
  meta = vcat([b.metadata for b in values(bundles)]...)
  meta = select(meta, [:Animal_nr, :Sex, :Genotype])
  filter!(row -> row.Animal_nr != 0, meta)
  rename!(meta, :Animal_nr => :Animal)
  meta.Group = string.(meta.Sex, "_", meta.Genotype)
  return meta
end

####################################################################################################

function load_experiments(params::DPParams)
  return load_experiments(params.expmeta, params.batches)
end

function load_experiments(metadata_path::String, batches::Vector{String})
  bundles = Dict{String,ExperimentBundle}()

  # TODO: add a branch to process csv files

  XLSX.openxlsx(metadata_path) do xf
    for sheetname in XLSX.sheetnames(xf)
      m = match(r"\d{4}-\d{2}-\d{2}", sheetname)
      isnothing(m) && continue
      datekey = m.match

      sheet = xf[sheetname]
      metadata = XLSX.gettable(sheet) |> DataFrame

      rename!(metadata, Symbol.(names(metadata)))
      metadata[!, :Cage_nr] = Int.(metadata[!, :Cage_nr])
      metadata[!, :Animal_nr] =
        map(x -> x == "EMPTY" ? 0 : parse(Int, string(x)), metadata[!, :Animal_nr])
      metadata[!, :Sex] = map(x -> x == "EMPTY" ? "" : string(x), metadata[!, :Sex])
      metadata[!, :Genotype] =
        map(x -> x == "EMPTY" ? "" : string(x), metadata[!, :Genotype])

      csvmatch = filter(f -> occursin(datekey, f), batches)
      if isempty(csvmatch)
        @warn "No CSV file found for date $datekey"
        continue
      end

      df = readdf(first(csvmatch); sep = ',')

      vars = setdiff(names(df), TIME_COLS)

      for c in TIME_COLS
        df[!, c] = DateTime.(df[!, c], dateformat"mm/dd/yyyy HH:MM:SS")
      end

      for c in vars
        df[!, c] = Float64.(df[!, c])
      end

      bundles[datekey] = ExperimentBundle(metadata, df)
    end
  end

  return bundles
end

####################################################################################################

function split_by_subject(bundle::ExperimentBundle; timecol::Symbol = Symbol(TIME_COLS[1]))
  suffix_ids = Dict{Symbol,Int}()
  for col in names(bundle.experiment)
    m = match(r"_(\d+)$", String(col))
    if m !== nothing
      suffix_ids[Symbol(col)] = parse(Int, m.captures[1])
    end
  end

  grouped = Dict{Int,Vector{Symbol}}()
  for (col, id) in suffix_ids
    push!(get!(grouped, id, Symbol[]), col)
  end

  raw_ids = bundle.metadata[!, ANIMAL_NR]
  idmap = Dict(i => raw_ids[i] for i = 1:nrow(bundle.metadata))

  dfs = Dict{Int,DataFrame}()
  for (suffix_id, cols) in grouped
    animal_nr = idmap[suffix_id]
    if animal_nr == 0
      continue
    end

    df = bundle.experiment[:, vcat([timecol], sort(cols))]
    renames = Dict(c => Symbol(replace(String(c), r"_\d+$" => "")) for c in cols)
    rename!(df, renames)
    bases = sort(Symbol.(replace.(String.(cols), r"_\d+$" => "")))
    dfs[animal_nr] = df[:, vcat([timecol], bases)]
  end

  return dfs
end

function split_by_subject(
  bundles::Dict{String,ExperimentBundle};
  timecol = Symbol(TIME_COLS[1]),
)
  consolidated = Dict{Int,DataFrame}()
  for bundle in values(bundles)
    dfs = split_by_subject(bundle; timecol = timecol)
    for (animal_nr, df) in dfs
      if haskey(consolidated, animal_nr)
        @warn "Duplicate Animal_nr $animal_nr across bundles; overwriting"
      end
      consolidated[animal_nr] = df
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
