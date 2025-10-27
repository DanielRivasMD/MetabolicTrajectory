####################################################################################################

using DataFrames
using Parameters: @with_kw

####################################################################################################

"""
    ExperimentBundle(metadata::DataFrame, experiment::DataFrame)

Container for a single experiment, holding:
- `metadata`  : DataFrame describing cages, animals, sex, and group (from XLSX sheet).
- `experiment`: DataFrame with the corresponding experimental time‑series data (from CSV).
"""
struct ExperimentBundle
  metadata::DataFrame
  experiment::DataFrame
end

####################################################################################################

"""
    TrajectoryParams(; kwargs...)

Configuration container for metabolic trajectory sampling and preprocessing.

# Fields
- `metadata::String`  
  Path to the XLSX file containing experiment metadata (cages, animals, sex, group).  
  **Required, no default.**

- `batches::Vector{String}`  
  One or more input files (CSV/other) containing metabolic data.  
  **Required, no default.**

- `nsamples::Int = 100`  
  Number of subsamples to draw per trajectory.

- `subsample_len::Float64 = 0.05`  
  Desired subsample length.  
  - If `0 < subsample_len < 1`, it is interpreted as a fraction of the full trajectory length (e.g. `0.05` = 5%).  
  - If `subsample_len >= 1`, it is interpreted as an absolute number of indices.

- `subsample_var::Float64 = 0.02`  
  Relative variance allowed in subsample length (e.g. +-2%).
"""
@with_kw mutable struct TrajectoryParams
  metadata::String
  batches::Vector{String}
  nsamples::Int = 100
  subsample_len::Float64 = 0.05
  subsample_var::Float64 = 0.01
end

####################################################################################################
