####################################################################################################

using Parameters: @with_kw

####################################################################################################

@with_kw mutable struct TrajectoryParams
  metadata::String
  batches::Vector{String}
  nsamples::Int = 100
  len::Float64 = 0.075
  var::Float64 = 0.01
end

####################################################################################################

function loadTparams(::Nothing, params::TrajectoryParams)
  return params
end

function loadTparams(path::String)
  params = TrajectoryParams()
  return loadTparams(path, params)
end

function loadTparams(path::String, params::TrajectoryParams)
  trajectory_cfg = path != "" ? symbolise_keys(TOML.parsefile(path)["trajectory"]) : Dict()
  return TrajectoryParams(; merge(struct_to_dict(params), trajectory_cfg)...)
end

####################################################################################################
