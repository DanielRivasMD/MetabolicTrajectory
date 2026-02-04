####################################################################################################
# cli args
####################################################################################################

# Load path definitions
include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

include(joinpath(Paths.UTIL, "args.jl"))      # Args API

# Parse CLI arguments
args = subsample_args()

####################################################################################################

include(joinpath(Paths.CONFIG, "dependencies.jl"))
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.CONFIG, "tparams.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "sigma.jl"))

####################################################################################################

include(joinpath(Paths.PIPE, "sigma.jl"))

####################################################################################################



####################################################################################################
