####################################################################################################

module DPREPL

####################################################################################################

using Avicenna.Flow: Cache, launch
using ..DPFlow: flow
using ..DPCore

####################################################################################################

export run

####################################################################################################

"""
    run(config_path::String; no_cache=false) -> Result

Run the data processing workflow from the REPL
"""
function run(config_path::String; no_cache = false)
  config = Dict{String,Any}("config_path" => config_path)
  cache = Cache("cache/data_processing", !no_cache)
  return launch(flow, config; cache = cache)
end

####################################################################################################

end

####################################################################################################
