####################################################################################################

module DPFlow

####################################################################################################

using Avicenna.Flow: Stage, Config
using ..DPCore

####################################################################################################

export flow

####################################################################################################

function ingest_stage(config::Dict, prev::Dict)
  config_path = config["config_path"]
  params = DPCore.load_ingest_params(config_path)
  # Allow overriding output directory from the upper config
  if haskey(config, "out_dir")
    params = DPCore.IngestParams(params.metadata, params.batches, config["out_dir"])
  end
  return DPCore.run_ingestion(params)
end

####################################################################################################

const flow = Config("data_processing", [Stage("01_ingest", ingest_stage, "1.0")], "1.0")

####################################################################################################

end

####################################################################################################
