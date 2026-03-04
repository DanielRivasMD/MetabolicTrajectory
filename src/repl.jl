####################################################################################################
# src/repl.jl
####################################################################################################

module REPL

####################################################################################################

using Avicenna.Workflow
using ..Process
using ..Workflows

####################################################################################################

function run_sigma(config_path::String; no_cache = false, plot = false)
  config = Core.load_toml_config(config_path)
  cache = Workflow.Cache("cache/sigma", !no_cache)
  workflow = plot ? Workflows.sigma_workflow_with_plots : Workflows.sigma_workflow
  return Workflow.run(workflow, config, cache = cache)
end

function plot_costmatrix(cost_path, ids_path, meta_path; kwargs...)
  # similar to CLI.plot_matrix but returns plot object
  cost = readdlm(cost_path, ',', Float64)
  ids_df = Core.readdf(ids_path; sep = ',')
  ids = Core.dataframe_to_ids(ids_df)
  meta = Core.readdf(meta_path; sep = ',')
  return Core.plot_grouped_costmatrix(cost, ids, meta; kwargs...)
end

####################################################################################################

end

####################################################################################################
