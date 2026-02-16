####################################################################################################

# Make src/pipe visible
push!(LOAD_PATH, @__DIR__)
using Pipeline

####################################################################################################

# Load the stage
include(joinpath(@__DIR__, "..", "bin", "mock.jl"))
using .Mock

####################################################################################################

# Create proper typed params
params = Mock.PipeParams()   # uses defaults

# Or load from TOML:
# params = Mock.loadPipeParams("config/pipeline.toml", Dict())

####################################################################################################

stages = [Stage("mock", Mock.run)]

run_pipeline(stages; params = params)

####################################################################################################
