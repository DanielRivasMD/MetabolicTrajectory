####################################################################################################

module Mock

####################################################################################################

export run, PipeParams, loadPipeParams

####################################################################################################

include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.UTIL, "args.jl"))
include(joinpath(Paths.CONFIG, "dependencies.jl"))
include(joinpath(Paths.UTIL, "pparams.jl"))

####################################################################################################

function run(params::PipeParams)
  @info "[mock] running with params: $(params)"

  if params.verbose
    @info "[mock] verbose mode enabled"
  end

  @info "[mock] processing input: $(params.input)"
  sleep(0.3)

  @info "[mock] writing output to: $(params.output)"
  sleep(0.3)

  @info "[mock] done"
end

end

####################################################################################################
# cli entry point
####################################################################################################

if abspath(PROGRAM_FILE) == @__FILE__
  args = Mock.pipe_args()
  params = Mock.loadPipeParams(get(args, "config", nothing), args)
  Mock.run(params)
end

####################################################################################################
