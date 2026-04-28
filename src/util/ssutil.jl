####################################################################################################

module SSCore

####################################################################################################

using DataFrames
using Dates
using Random
using ..USCore: readdf, writedf

####################################################################################################

struct SubSamplingParams
  input_dir::String
  nsamples::Int
  len::Float64            # relative if < 1, absolute if >= 1
  var::Float64            # relative jitter (proportion of length)
  limits::Tuple{Int,Int}  # (lo, hi) 1‑based indices; (0,0) -> no limit
  seed::Union{Int,Nothing}
  variables::Union{Nothing,Vector{String}}
end

####################################################################################################

function load_existing_data(input_dir::String)
  meta_path = joinpath(input_dir, "meta.csv")
  isfile(meta_path) || error("meta.csv not found in $input_dir")
  meta = readdf(meta_path; sep = ',')

  animals = Dict{Int,DataFrame}()
  for file in readdir(input_dir)
    m = match(r"exp_(\d+)\.csv", file)
    m === nothing && continue
    id = parse(Int, m.captures[1])
    df = readdf(joinpath(input_dir, file); sep = ',')
    animals[id] = df
  end
  return meta, animals
end

####################################################################################################

function generate_subsamples(
  signal::Vector{Float64},
  nsamples::Int,
  len::Float64,
  var::Float64,
  limits::Tuple{Int,Int};
  seed::Union{Int,Nothing} = nothing,
)
  lo, hi = limits
  if lo != 0 || hi != 0
    lo = max(1, lo)
    hi = min(hi, length(signal))
    signal = signal[lo:hi]
  else
    lo = 1
    hi = length(signal)
  end

  n = length(signal)
  if n == 0
    return Tuple{Vector{Float64},Int,Int}[]
  end

  base_len = len < 1.0 ? round(Int, len * n) : round(Int, len)

  rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)

  results = Vector{Tuple{Vector{Float64},Int,Int}}()
  for _ = 1:nsamples
    jitter = round(Int, base_len * (1 + (rand(rng) * var * 2 - var)))
    actual_len = clamp(jitter, 1, n)
    if n <= actual_len
      continue
    end
    s = rand(rng, 1:(n-actual_len))
    e = s + actual_len - 1
    push!(results, (signal[s:e], lo + s - 1, lo + e - 1))
  end
  return results
end

####################################################################################################

subsample_id(animal, start_idx) = animal * 100_000 + start_idx

####################################################################################################

function write_subsample_csv(
  out_dir::String,
  unique_id::Int,
  signal::Vector{Float64},
  var_name::String,
)
  df = DataFrame(
    Date_Time = 0:length(signal)-1,  # simple index
    (Symbol(var_name) => signal),
  )
  writedf(joinpath(out_dir, "exp_$(unique_id).csv"), df; sep = ',')
end

####################################################################################################

function run_subsampling(params::SubSamplingParams)
  in_dir = params.input_dir
  meta, animals = load_existing_data(in_dir)

  # Determine variable list
  first_animal = first(values(animals))
  all_vars = setdiff(names(first_animal), ["Date_Time"])
  vars_to_process =
    isnothing(params.variables) ? all_vars : intersect(params.variables, all_vars)

  # Output directory
  parent = dirname(in_dir)
  base = basename(in_dir)
  seed_tag = isnothing(params.seed) ? "XX" : string(params.seed)
  out_dir = joinpath(parent, "sub_$(base)_$(seed_tag)")
  mkpath(out_dir)

  # New metadata DataFrame
  new_meta = DataFrame(
    subsample_id = Int[],
    original_animal = Int[],
    start_idx = Int[],
    end_idx = Int[],
    variable = String[],
    Sex = String[],
    Genotype = String[],
    Group = String[],
  )

  # Lookup for metadata of original animals
  meta_dict =
    Dict(row.Animal => (row.Sex, row.Genotype, row.Group) for row in eachrow(meta))

  for var in vars_to_process
    for (animal_id, df) in animals
      signal = Float64.(collect(skipmissing(df[!, var])))
      if isempty(signal)
        continue
      end

      subs = generate_subsamples(
        signal,
        params.nsamples,
        params.len,
        params.var,
        params.limits;
        seed = params.seed,
      )

      for (sub_signal, start_idx, end_idx) in subs
        uid = subsample_id(animal_id, start_idx)
        write_subsample_csv(out_dir, uid, sub_signal, var)
        sex, geno, group = meta_dict[animal_id]
        push!(new_meta, (uid, animal_id, start_idx, end_idx, var, sex, geno, group))
      end
    end
  end

  # Save new metadata
  writedf(joinpath(out_dir, "meta_subsamples.csv"), new_meta; sep = ',')

  return Dict("output_dir" => out_dir, "n_subsamples" => nrow(new_meta))
end

####################################################################################################

function load_toml_config(path::String)
  raw = TOML.parsefile(path)
  return raw["subsampling"]
end

####################################################################################################

end

####################################################################################################
