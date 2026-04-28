####################################################################################################

module HCREPL

####################################################################################################

using Avicenna.Flow: Cache, launch
using ..HCFlow: flow
using ..HCCore

####################################################################################################

export run

####################################################################################################

function run(
  matrix_path::String,
  meta_path::String;
  output_dir::String = "clustering",
  k::Int = 2,
  title::String = "",
  no_cache::Bool = false,
)
  config = Dict{String,Any}(
    "matrix" => matrix_path,
    "meta" => meta_path,
    "output_dir" => output_dir,
    "k" => k,
    "title" => title,
  )
  cache = Cache("cache/hierarchical_clustering", !no_cache)
  return launch(flow, config; cache = cache)
end

####################################################################################################

end

####################################################################################################
