####################################################################################################

using Parameters: @with_kw
using TOML

####################################################################################################

# CLI > TOML > defaults
function loadParams(
  path::Union{Nothing,String},
  ::Type{T},
  args = nothing;
  section = nothing,
  postprocess = nothing,
) where {T}

  base = struct_to_dict(T())

  if !(path === nothing || path == "")
    raw = TOML.parsefile(path)
    cfg = section === nothing ? raw : raw[string(section)]
    cfg = symbolise_keys(cfg)

    if postprocess !== nothing
      cfg = postprocess(cfg)
    end

    base = merge(base, cfg)
  end

  if args !== nothing
    valid = Set(keys(base))
    cli_cfg = Dict(Symbol(k) => v for (k, v) in args if Symbol(k) in valid && v !== nothing)
    base = merge(base, cli_cfg)
  end

  return T(; base...)
end

####################################################################################################

@with_kw mutable struct TrajectoryParams
  metadata::String = ""
  batches::Vector{String} = []
  out::String = ""
  nsamples::Int = 10
  len::Float64 = 0.075
  var::Float64 = 0.005
  limits::Tuple{Int,Int} = (0, 0)
end

####################################################################################################

function loadTrajectoryParams(path::Union{Nothing,String}, args::Dict)
  return loadParams(
    path,
    TrajectoryParams,
    args;
    section = :trajectory,
    postprocess = cfg -> begin
      if :limit1 in keys(cfg) && :limit2 in keys(cfg)
        cfg[:limits] = (cfg[:limit1], cfg[:limit2])
        delete!(cfg, :limit1)
        delete!(cfg, :limit2)
      end
      cfg
    end,
  )
end

####################################################################################################

@with_kw mutable struct PipeParams
  input::String = "iparmas"
  output::String = "oparams"
  verbose::Bool = false
end

####################################################################################################

function loadPipeParams(config_path::Union{Nothing,String}, args)
  return loadParams(config_path, PipeParams, args; section = :pipe)
end

####################################################################################################
