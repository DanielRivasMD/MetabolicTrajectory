####################################################################################################

using Parameters: @with_kw
using TOML

####################################################################################################

"""
    PipeParams(; kwargs...)

Parameter container for the mock stage.

# Fields
- `input::String`
- `output::String`
- `verbose::Bool = false`
"""
@with_kw mutable struct PipeParams
  input::String = "iparmas"
  output::String = "oparams"
  verbose::Bool = false
end

####################################################################################################

"""
    loadPipeParams(config_path::Union{Nothing,String}, args)

Load defaults from PipeParams, optionally override from a TOML file,
and finally override with CLI args.
"""
function loadPipeParams(config_path::Union{Nothing,String}, args)
  params = PipeParams()

  toml_cfg = Dict()
  if config_path !== nothing
    raw = TOML.parsefile(config_path)
    if haskey(raw, "pipe")
      toml_cfg = symbolise_keys(raw["pipe"])
    end
  end

  merged = merge(struct_to_dict(params), toml_cfg)
  valid = Set(propertynames(params))
  cli_cfg = Dict(Symbol(k) => v for (k, v) in args if Symbol(k) in valid && v !== nothing)
  merged = merge(merged, cli_cfg)

  return PipeParams(; merged...)
end

####################################################################################################
