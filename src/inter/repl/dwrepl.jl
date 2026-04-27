####################################################################################################

module DWREPL

####################################################################################################

using Avicenna.Flow: Cache, launch
using ..DWFlow: flow
using ..DWCore

####################################################################################################

export run

####################################################################################################

"""
    run(input_dir::String; output_dir="dtw_full", variables=nothing, no_cache=false) -> Result

Run the DTW analysis workflow from the REPL.
"""
function run(
  input_dir::String;
  output_dir::String = "dtw_full",
  variables::Union{Nothing,Vector{String}} = nothing,
  no_cache::Bool = false,
)
  config = Dict{String,Any}("input_dir" => input_dir, "output_dir" => output_dir)
  if variables !== nothing
    config["variables"] = variables
  end
  cache = Cache("cache/dynamic_time_warping", !no_cache)
  return launch(flow, config; cache = cache)
end

####################################################################################################

end

####################################################################################################
