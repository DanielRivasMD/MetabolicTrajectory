####################################################################################################

push!(LOAD_PATH, @__DIR__)
using Pipeline

####################################################################################################

include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.BIN, "sigma.jl"));
using .Sigma;
include(joinpath(Paths.BIN, "saturated.jl"));
using .Saturated;

####################################################################################################

sigma_params = Sigma.TrajectoryParams()

####################################################################################################

stages = [
  Stage("sigma", Sigma.run)
  Stage("saturated", Saturated.run)
]

run_pipeline(stages; params = Dict())

####################################################################################################
