####################################################################################################

module Paths

"project"
const PROJECT = normpath(joinpath(@__DIR__, "..", ".."))

# tmp csv
const TMP = joinpath(PROJECT, "tmp")

"data"
const DATA = joinpath(PROJECT, "data")
const SIGMA = joinpath(DATA, "sigma")
const SCSV = joinpath(SIGMA, "csv")
const STMP = joinpath(SIGMA, "tmp")
const SXLSX = joinpath(SIGMA, "xlsx")
const HMGCR = joinpath(DATA, "HMGCR_TAM")
const HCSV = joinpath(HMGCR, "csv")
const HXLSX = joinpath(HMGCR, "xlsx")

"graph"
const GRAPH = joinpath(PROJECT, "graph")
const HTML = joinpath(GRAPH, "html")

"src"
const SRC = joinpath(PROJECT, "src")
const BIN = joinpath(SRC, "bin")
const CONFIG = joinpath(SRC, "config")
const LIVE = joinpath(SRC, "live")
const UTIL = joinpath(SRC, "util")

"Ensure directories exist (for outputs)"
function ensure_dirs()
  for d in (DATA, SIGMA, SCSV, SXLSX, HMGCR, HCSV, HXLSX, GRAPH, HTML)
    isdir(d) || mkpath(d)
  end
end

end

####################################################################################################
