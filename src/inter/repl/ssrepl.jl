# src/inter/repl/ssrepl.jl
module SSREPL

using Avicenna.Flow: Cache, launch
using ..SSFlow: flow
using ..SSCore

export run

function run(config_path::String; no_cache::Bool = false)
  flow_config = Dict{String,Any}("config_path" => config_path)
  cache = Cache("cache/ss", !no_cache)
  return launch(flow, flow_config; cache = cache)
end

end
