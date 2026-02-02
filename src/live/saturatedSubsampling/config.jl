####################################################################################################

using DataFrames
using Dates
using Parameters: @with_kw

####################################################################################################

struct ExperimentBundle
  metadata::DataFrame
  experiment::DataFrame
end

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

@with_kw mutable struct TrajectoryParams
  metadata::String
  batches::Vector{String}
  nsamples::Int = 100
  subsample_len::Float64 = 0.075
  subsample_var::Float64 = 0.0
end

####################################################################################################
