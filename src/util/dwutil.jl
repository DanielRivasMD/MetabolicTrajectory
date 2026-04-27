####################################################################################################

module DWCore

####################################################################################################

using DataFrames
using DelimitedFiles
using DynamicAxisWarping
using Distances
using Statistics
using ..USCore: readdf, writedf

####################################################################################################

const TIME_COL = "Date_Time"

####################################################################################################

"""
    load_experiment_data(input_dir::String) -> (meta::DataFrame, animals::Dict{Int,DataFrame}, variables::Vector{String})

Loads per‑animal CSV files and meta.csv from `input_dir`.
Returns meta DataFrame, dictionary mapping animal ID to its time‑series DataFrame,
and a list of measurement variable names (all columns except `TIME_COL`).
"""
function load_experiment_data(input_dir::String)
  # Load metadata
  meta_path = joinpath(input_dir, "meta.csv")
  if !isfile(meta_path)
    error("meta.csv not found in $input_dir")
  end
  meta = readdf(meta_path; sep = ',')

  # Load per‑animal files
  animals = Dict{Int,DataFrame}()
  for file in readdir(input_dir)
    m = match(r"exp_(\d+)\.csv", file)
    m === nothing && continue
    id = parse(Int, m.captures[1])
    df = readdf(joinpath(input_dir, file); sep = ',')
    animals[id] = df
  end

  # Extract variable names from first animal
  first_df = first(values(animals))
  variables = setdiff(names(first_df), [TIME_COL])

  return meta, animals, variables
end

####################################################################################################

"""
    dtw_matrix_for_variable(animals::Dict{Int,DataFrame}, var::String) -> (subjects::Vector{Int}, matrix::Matrix{Float64})

Computes a normalised DTW cost matrix between all animals for the given variable.
The matrix is indexed by the sorted animal IDs.
"""
function dtw_matrix_for_variable(animals::Dict{Int,DataFrame}, var::String)
  subjects = sort(collect(keys(animals)))
  N = length(subjects)

  # Extract signals
  signals = Dict{Int,Vector{Float64}}()
  for subj in subjects
    df = animals[subj]
    if var in names(df)
      sig = collect(skipmissing(df[!, var]))
      if !isempty(sig)
        signals[subj] = sig
      end
    end
  end

  # Keep only subjects that have data
  valid = sort(collect(keys(signals)))
  Nvalid = length(valid)
  matrix = zeros(Float64, Nvalid, Nvalid)

  for i = 1:Nvalid, j = i:Nvalid
    s1 = signals[valid[i]]
    s2 = signals[valid[j]]
    cost, _, _ = dtw(s1, s2, SqEuclidean())
    norm_cost = cost / mean([length(s1), length(s2)])
    matrix[i, j] = norm_cost
    matrix[j, i] = norm_cost
  end

  return valid, matrix
end

####################################################################################################

"""
    run_dtw(input_dir::String, output_dir::String; variables::Union{Nothing,Vector{String}}=nothing) -> Dict

Performs full DTW analysis:
  1. Load data from `input_dir`.
  2. For each variable (or only those in `variables`):
       - compute normalised DTW matrix
       - save matrix CSV and order IDs CSV to `output_dir`.
Returns a summary dictionary.
"""
function run_dtw(
  input_dir::String,
  output_dir::String;
  variables::Union{Nothing,Vector{String}} = nothing,
)
  meta, animals, all_vars = load_experiment_data(input_dir)
  mkpath(output_dir)

  vars_to_process = isnothing(variables) ? all_vars : intersect(variables, all_vars)
  results = Dict{String,Dict}()

  for var in vars_to_process
    subjects, matrix = dtw_matrix_for_variable(animals, var)

    # Save matrix with subject IDs as headers
    out_matrix_path = joinpath(output_dir, "$(var)_dtw_matrix.csv")
    save_dtw_matrix(out_matrix_path, subjects, matrix)

    # Save subject order list (just in case)
    out_ids_path = joinpath(output_dir, "$(var)_subjects.csv")
    df_ids = DataFrame(subject = subjects)
    writedf(out_ids_path, df_ids; sep = ',')

    results[var] = Dict("n_subjects" => length(subjects), "matrix_file" => out_matrix_path)
  end

  return Dict(
    "output_dir" => output_dir,
    "variables_processed" => collect(keys(results)),
    "details" => results,
  )
end

# Helper: write matrix with subject heading
function save_dtw_matrix(path::String, subjects::Vector{Int}, matrix::Matrix{Float64})
  N = length(subjects)
  out = Array{Any}(undef, N + 1, N + 1)
  out[1, 1] = "subject"
  for i = 1:N
    out[1, i+1] = subjects[i]
    out[i+1, 1] = subjects[i]
    for j = 1:N
      out[i+1, j+1] = matrix[i, j]
    end
  end
  writedlm(path, out, ',')
  return nothing
end

####################################################################################################

end

####################################################################################################
