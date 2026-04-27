####################################################################################################

module HCCLI

####################################################################################################

using ArgParse
using Avicenna.Flow: Cache, launch
using ..HCFlow: flow

####################################################################################################

export run

####################################################################################################

function run(args::Vector{String})
  s = ArgParseSettings()
  @add_arg_table! s begin
    "--matrix", "-m"
    help = "Path to DTW cost matrix CSV (with subject header)"
    required = true
    "--meta"
    help = "Path to meta CSV"
    required = true
    "--output-dir", "-o"
    help = "Output directory for plots and assignments"
    default = "clustering"
    "--k"
    help = "Number of clusters (default: 4)"
    arg_type = Int
    default = 4
    "--title"
    help = "Plot title"
    default = ""
    "--no-cache"
    help = "Disable caching"
    action = :store_true
  end
  parsed = parse_args(args, s)

  config = Dict{String,Any}(
    "matrix" => parsed["matrix"],
    "meta" => parsed["meta"],
    "output_dir" => parsed["output-dir"],
    "k" => parsed["k"],
    "title" => parsed["title"],
  )

  cache = Cache("cache/hc", !parsed["no-cache"])
  result = launch(flow, config; cache = cache)
  summary = result.stage_outputs["02_cluster"]
  println("Clustering complete. Results in: ", summary["output_dir"])
  return 0
end

####################################################################################################

end

####################################################################################################
