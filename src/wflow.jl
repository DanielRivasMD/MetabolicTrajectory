####################################################################################################
# src/workflows.jl
####################################################################################################

module Workflows

####################################################################################################

using Avicenna.Workflow
using ..Process

####################################################################################################
# Stage 1: Load experiments and split by animal
####################################################################################################

function load_and_split_stage(config, _)
  params = Process.TrajectoryParams(; config...)   # convert Dict to struct
  bundles = Process.load_experiments(params)       # (you need to implement this)
  subdfs = Process.split_by_animal(bundles)
  # Save per‑animal CSVs (optional – maybe you want to cache this)
  Process.writedf_dict(params.meta_path, subdfs)
  return subdfs
end

####################################################################################################
# Stage 2: For each variable, collect subsamples and compute DTW
####################################################################################################

function dtw_stage(config, prev)
  subdfs = prev["load_and_split"]
  params = Process.TrajectoryParams(; config...)

  first_df = first(values(subdfs))
  vars = Symbol.(names(first_df)[2:end])

  results = Dict{Symbol,Any}()
  for var in vars
    container = Process.collect_subsamples(subdfs, var, params, true)  # diff=true
    isempty(container.subsamples) && continue

    cost_matrix = Process.compute_dtw_matrix(container.subsamples)

    # Save to disk (or you could store in the result for later)
    tag = Process.matrix_tag(params)
    outdir = joinpath(params.dtw_path, tag)
    mkpath(outdir)
    Process.writedf(joinpath(outdir, "$(var)_cost_matrix.csv"), cost_matrix)
    Process.writedf(
      joinpath(outdir, "$(var)_order_ids.csv"),
      Process.ids_to_dataframe(container.ids),
    )

    results[var] = (matrix = cost_matrix, ids = container.ids)
  end
  return results
end

####################################################################################################
# Stage 3: Generate plots (optional)
####################################################################################################

function plot_stage(config, prev)
  # Use prev["dtw"] to create plots, save them, and return plot objects
  # ...
  return Dict("plot_paths" => ["..."])
end

####################################################################################################
# Assemble workflows
####################################################################################################

const sigma_workflow = Workflow.WorkflowConfig(
  "sigma_analysis",
  [
    Workflow.Stage("load_and_split", load_and_split_stage, "1.0"),
    Workflow.Stage("dtw", dtw_stage, "1.0"),
  ],
  "1.0",
)

const sigma_workflow_with_plots = Workflow.WorkflowConfig(
  "sigma_analysis_plots",
  [
    Workflow.Stage("load_and_split", load_and_split_stage, "1.0"),
    Workflow.Stage("dtw", dtw_stage, "1.0"),
    Workflow.Stage("plot", plot_stage, "1.0"),
  ],
  "1.0",
)

####################################################################################################

end

####################################################################################################
