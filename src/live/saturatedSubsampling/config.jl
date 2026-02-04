####################################################################################################

using Parameters: @with_kw

####################################################################################################

@with_kw mutable struct TrajectoryParams
  metadata::String
  batches::Vector{String}
  nsamples::Int = 100
  subsample_len::Float64 = 0.030
  subsample_var::Float64 = 0.005
  sample::Tuple{Int, Int} = (0, 0)
end

####################################################################################################

using Parameters: @with_kw

####################################################################################################

using TOML
using FilePathsBase: basename, splitext, joinpath, dirname, isabspath, isdir, mkpath

####################################################################################################

function loadTparams(path::String)
  params = TrajectoryParams()
  sample_cfg = path != "" ? symbolise_keys(TOML.parsefile(path)["sample"]) : Dict()
  return TrajectoryParams(; merge(struct_to_dict(params), sample_cfg)...)
end

####################################################################################################
