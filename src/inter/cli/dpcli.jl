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
    "--config"
    help = "TOML configuration file"
    required = true
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
  cache = Cache("cache/data_processing", !parsed["no-cache"])
  result = launch(flow, flow_config; cache = cache)

  println("Data processing completed")
  params = result.stage_outputs["01_load_params"]
  println("Output directory: ", params.outdir)
  return 0
end

####################################################################################################

end

####################################################################################################
