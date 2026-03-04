#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "MetabolicTrajectories.jl"))
using .MetabolicTrajectories

MetabolicTrajectories.CLI.main(ARGS)
