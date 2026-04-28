####################################################################################################

module HCFlow

####################################################################################################

using Avicenna.Flow: Stage, Config
using ..HCCore

####################################################################################################

export flow

####################################################################################################

function load_config_stage(config::Dict)
  return Dict(
    "matrix_path" => config["matrix"],
    "meta_path" => config["meta"],
    "output_dir" => config["output_dir"],
    "k" => get(config, "k", 2),
    "title" => get(config, "title", ""),
  )
end

####################################################################################################

function cluster_stage(prev::Dict)
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
