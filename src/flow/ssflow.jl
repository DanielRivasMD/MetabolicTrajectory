# src/flow/ssflow.jl
module SSFlow

using Avicenna.Flow: Stage, Config
using ..SSCore

export flow

function load_config_stage(config::Dict, prev::Dict)
  # config should contain "config_path" pointing to TOML
  toml = SSCore.load_toml_config(config["config_path"])   # to be defined in SSCore
  return Dict(
    "params" => SSCore.SubSamplingParams(
      toml["input_dir"],
      toml["nsamples"],
      toml["len"],
      toml["var"],
      get(toml, "limits", (0, 0)),
      get(toml, "seed", nothing),
      get(toml, "variables", nothing),
    ),
  )
end

function run_subsampling_stage(config::Dict, prev::Dict)
  params = prev["01_load_config"]["params"]
  return SSCore.run_subsampling(params)
end

# Add the TOML loading helper in SSCore
# We'll update SSCore to include load_toml_config

const flow = Config(
  "subsampling",
  [
    Stage("01_load_config", load_config_stage, "1.0"),
    Stage("02_subsample", run_subsampling_stage, "1.0"),
  ],
  "1.0",
)

end
