####################################################################################################

module DPFlow

####################################################################################################

using Avicenna.Flow: Stage, Config
using ..DPCore

####################################################################################################

export flow

####################################################################################################

const flow = Config(
  "data_processing",
  [
    Stage(
      "01_load_params",
      (config, _) -> begin
        return DPCore.load_params(config["config_path"])
      end,
      "1.0",
    ),
    Stage("02_make_path", (_, prev) -> begin
      params = prev["01_load_params"]
      mkpath(params.outdir)
      return nothing
    end, "1.0"),
    Stage(
      "03_load_experiments",
      (_, prev) -> begin
        params = prev["01_load_params"]
        return DPCore.load_experiments(params)
      end,
      "1.0",
    ),
    Stage("04_build_meta", (_, prev) -> begin
      bundles = prev["03_load_experiments"]
      return DPCore.build_meta(bundles)
    end, "1.0"),
    Stage(
      "05_write_meta",
      (_, prev) -> begin
        params = prev["01_load_params"]
        meta = prev["04_build_meta"]
        DPCore.writedf(joinpath(params.outdir, params.outmeta), meta; sep = ',')
        return nothing

      end,
      "1.0",
    ),
    Stage(
      "06_write_subjects",
      (_, prev) -> begin
        params = prev["01_load_params"]
        bundles = prev["03_load_experiments"]
        dfs = DPCore.split_by_subject(bundles)
        DPCore.writedf_dict(params.outdir, dfs; sep = ',')
        return nothing
      end,
      "1.0",
    ),
  ],
  "1.0",
)

####################################################################################################

end

####################################################################################################
