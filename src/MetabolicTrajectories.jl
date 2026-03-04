# src/MetabolicTrajectories.jl
module MetabolicTrajectories

# Package for trajectory analysis using Avicenna framework

# Export symbols from submodules for convenience
export Process, Workflows, CLI, REPL, Document
export run_sigma, plot_costmatrix  # example re-exports from REPL/CLI

# -----------------------------------------------------------------------------
# Include utility modules first (they are used by others)
# -----------------------------------------------------------------------------
include("util/ioDataFrame.jl")
include("util/params.jl")
include("util/wrangle.jl")

# -----------------------------------------------------------------------------
# Core module – pure analysis logic
# -----------------------------------------------------------------------------
include("core.jl")

# -----------------------------------------------------------------------------
# Workflow definitions using Avicenna.Workflow
# -----------------------------------------------------------------------------
include("wflow.jl")

# -----------------------------------------------------------------------------
# Interfaces: REPL, CLI, Document
# -----------------------------------------------------------------------------
include("repl.jl")
include("cli.jl")
include("document.jl")

# -----------------------------------------------------------------------------
# Re-export commonly used functions (optional)
# -----------------------------------------------------------------------------
# From Core
export TrajectoryParams, SubSampleID, SubSampleContainer
export collect_subsamples, compute_dtw_matrix, diff_cumulative
export readdf, writedf, readdf_dict, writedf_dict  # from utils/io

# From Workflows
export sigma_workflow, sigma_workflow_with_plots

# From REPL
export run_sigma, plot_costmatrix

# From CLI – typically not re-exported, but you could
# export main if you want to call it programmatically

end # module
