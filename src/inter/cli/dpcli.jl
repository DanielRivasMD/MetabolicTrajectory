####################################################################################################

module DPCLI

####################################################################################################

using ArgParse
using Avicenna.Flow: Cache, launch
using ..DPFlow: flow

####################################################################################################

export run

####################################################################################################

function run(args::Vector{String})
  s = ArgParseSettings()
  @add_arg_table! s begin
    "--config", "-c"
    help = "Path to TOML configuration file for ingestion"
    required = true
    "--out-dir", "-o"
    help = "Override default output directory"
    default = nothing
    "--no-cache"
    help = "Disable caching"
    action = :store_true
  end
  parsed = parse_args(args, s)

  config_path = parsed["config"]
  if !isfile(config_path)
    @error "Configuration file not found: $config_path"
    return 1
  end

  flow_config = Dict{String,Any}("config_path" => config_path)
  if parsed["out-dir"] !== nothing
    flow_config["out_dir"] = parsed["out-dir"]
  end

  cache = Cache("cache/data_processing", !parsed["no-cache"])
  result = launch(flow, flow_config; cache = cache)

  println("Data processing completed.")
  println("Output directory: ", result.stage_outputs["01_ingest"]["out_dir"])
  return 0
end

####################################################################################################

end

####################################################################################################
