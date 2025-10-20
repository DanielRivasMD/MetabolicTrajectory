####################################################################################################

module Paths

"project"
const PROJECT = normpath(joinpath(@__DIR__, "..", ".."))

"data"
const DATA = joinpath(PROJECT, "data")
const CSV = joinpath(DATA, "csv")
const HMGCR = joinpath(DATA, "HMGCR_TAM")
const XLSX = joinpath(DATA, "xlsx")

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
  for d in (DATA, XLSX, GRAPH)
    isdir(d) || mkpath(d)
  end
end

end

####################################################################################################
