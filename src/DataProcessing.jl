####################################################################################################

"""
Performs data processing:
  1. Load experiments from metadata XLSX and CSV batches
  2. Build a unified metadata DataFrame (Animal, Sex, Genotype)
  3. Split time‑series data by animal
  4. Write meta.csv and per‑animal CSVs to `params.outdir`
"""
module DP

####################################################################################################

include("util/shared.jl")
include("util/dputil.jl")
include("flow/dpflow.jl")
include("inter/cli/dpcli.jl")
include("inter/repl/dprepl.jl")

####################################################################################################

export DPCore, DPFlow, DPCLI, DPREPL

####################################################################################################

end

####################################################################################################
