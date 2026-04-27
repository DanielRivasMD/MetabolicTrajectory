####################################################################################################

module HCFlow

####################################################################################################

using Avicenna.Flow: Stage, Config
using ..HCCore

####################################################################################################

export flow

####################################################################################################

"""
Stage 01: Load DTW matrix and metadata, returns paths and optional k.
"""
function load_config_stage(config::Dict, prev::Dict)
  return Dict(
    "matrix_path" => config["matrix"],
    "meta_path" => config["meta"],
    "output_dir" => config["output_dir"],
    "k" => get(config, "k", 4),
    "title" => get(config, "title", ""),
  )
end

####################################################################################################

"""
Stage 02: Perform clustering and produce outputs.
"""
function cluster_stage(config::Dict, prev::Dict)
  args = prev["01_load_config"]
  return HCCore.cluster_analysis(
    args["matrix_path"],
    args["meta_path"],
    args["output_dir"];
    k = args["k"],
    plot_title = args["title"],
  )
end

####################################################################################################

const flow = Config(
  "hierarchical_clustering",
  [
    Stage("01_load_config", load_config_stage, "1.0"),
    Stage("02_cluster", cluster_stage, "1.0"),
  ],
  "1.0",
)

####################################################################################################

end

####################################################################################################
