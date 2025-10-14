####################################################################################################

module Paths

"Project root"
const PROJECT = normpath(joinpath(@__DIR__, "..", ".."))

"Data"
const DATA = joinpath(PROJECT, "data")
const XLSX = joinpath(DATA, "xlsx")

"Graph"
const GRAPH = joinpath(PROJECT, "graph")

"Source"
const SRC = joinpath(PROJECT, "src")
const CONFIG = joinpath(SRC, "config")
const UTIL = joinpath(SRC, "util")
const BIN = joinpath(SRC, "bin")

"Ensure directories exist (for outputs)"
function ensure_dirs()
  for d in (DATA, XLSX, GRAPH)
    isdir(d) || mkpath(d)
  end
end

end

####################################################################################################
