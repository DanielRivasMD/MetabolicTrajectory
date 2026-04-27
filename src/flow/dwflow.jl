####################################################################################################

module DWFlow

####################################################################################################

using Avicenna.Flow: Stage, Config
using ..DWCore

####################################################################################################

export flow

####################################################################################################

"""
Stage 01: Load experiment data from input directory.
Config must contain `"input_dir"`.
"""
function load_data_stage(config::Dict, prev::Dict)
  input_dir = config["input_dir"]
  meta, animals, variables = DWCore.load_experiment_data(input_dir)
  return Dict("meta" => meta, "animals" => animals, "variables" => variables)
end

"""
Stage 02: Compute DTW matrices for selected variables and save to output directory.
Config must contain `"output_dir"` and optionally `"variables"` (list of variable names).
"""
function compute_dtw_stage(config::Dict, prev::Dict)
  data = prev["01_load_data"]
  animals = data["animals"]
  variables = data["variables"]
  output_dir = config["output_dir"]
  selected = get(config, "variables", nothing)
  return DWCore.run_dtw("", output_dir; variables = selected)  # we only need animals, so pass dummy input_dir
  # Actually run_dtw expects input_dir; we can create a wrapper that accepts already loaded data.
end

# A wrapper that accepts pre‑loaded data
function run_dtw_with_data(meta, animals, variables, output_dir; selected = nothing)
  vars_to_process = isnothing(selected) ? variables : intersect(selected, variables)
  mkpath(output_dir)
  results = Dict()
  for var in vars_to_process
    subjects, matrix = DWCore.dtw_matrix_for_variable(animals, var)
    out_matrix = joinpath(output_dir, "$(var)_dtw_matrix.csv")
    DWCore.save_dtw_matrix(out_matrix, subjects, matrix)
    results[var] = (n_subjects = length(subjects), matrix_file = out_matrix)
  end
  return Dict(
    "output_dir" => output_dir,
    "variables_processed" => collect(keys(results)),
    "details" => results,
  )
end

####################################################################################################

const flow = Config(
  "dynamic_time_warping",
  [
    Stage("01_load_data", load_data_stage, "1.0"),
    Stage(
      "02_compute_dtw",
      (config, prev) -> begin
        data = prev["01_load_data"]
        return run_dtw_with_data(
          data["meta"],
          data["animals"],
          data["variables"],
          config["output_dir"];
          selected = get(config, "variables", nothing),
        )
      end,
      "1.0",
    ),
  ],
  "1.0",
)

####################################################################################################

end

####################################################################################################
