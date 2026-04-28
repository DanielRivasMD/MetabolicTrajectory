####################################################################################################

module SSCLI

####################################################################################################

using ArgParse
using Avicenna.Flow: Cache, launch
using ..SSFlow: flow

####################################################################################################

export run

####################################################################################################

function run(args::Vector{String})
  s = ArgParseSettings()
  @add_arg_table! s begin
    "--config", "-c"
    help = "Path to TOML configuration file for subsampling"
    required = true
    "--no-cache"
    help = "Disable caching"
    action = :store_true
  end
  parsed = parse_args(args, s)

  config_path = parsed["config"]
  isfile(config_path) || error("Config file not found: $config_path")

  flow_config = Dict{String,Any}("config_path" => config_path)
  cache = Cache("cache/ss", !parsed["no-cache"])
  result = launch(flow, flow_config; cache = cache)

  summary = result.stage_outputs["02_subsample"]
  println("Subsampling finished. Output directory: ", summary["output_dir"])
  return 0
end

####################################################################################################

end

####################################################################################################
