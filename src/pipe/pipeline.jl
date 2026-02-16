# Make src/pipe visible
push!(LOAD_PATH, @__DIR__)

# Load the orchestrator
using Pipeline

# Load stage modules manually
include("stages/parse.jl")
include("stages/fft.jl")
include("stages/report.jl")

using .ParseStage
using .FFTStage
using .ReportStage

params = Dict("input" => "data/", "output" => "results/")

stages = [Stage("parse", run_parse), Stage("fft", run_fft), Stage("report", run_report)]

run_pipeline(stages; params = params)
