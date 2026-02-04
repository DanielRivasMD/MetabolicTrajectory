####################################################################################################

using Parameters: @with_kw

####################################################################################################

@with_kw mutable struct TrajectoryParams
  metadata::String
  batches::Vector{String}
  nsamples::Int = 10
  len::Float64 = 0.075
  var::Float64 = 0.005
  limits::Tuple{Int,Int} = (0, 0)
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

  # rebuild tuple
  if :limit1 in keys(trajectory_cfg) && :limit2 in keys(trajectory_cfg)
    trajectory_cfg[:limits] = (trajectory_cfg[:limit1], trajectory_cfg[:limit2])
    delete!(trajectory_cfg, :limit1)
    delete!(trajectory_cfg, :limit2)
  end

  return TrajectoryParams(; merge(struct_to_dict(params), trajectory_cfg)...)
end

####################################################################################################

function tparam_tag(p::TrajectoryParams)
  return @sprintf(
    "n%d_len%.3f_var%.3f_limits%d_%d",
    p.nsamples,
    p.len,
    p.var,
    p.limits[1],
    p.limits[2]
  )
end

####################################################################################################
