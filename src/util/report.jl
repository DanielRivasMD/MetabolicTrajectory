####################################################################################################

include(joinpath(@__DIR__, "..", "config", "paths.jl"))
using .Paths

####################################################################################################

using Weave
using DataFrames
using PrettyTables

####################################################################################################

function report(root::String)
  input_file = joinpath(Paths.SRC, "live", root * ".jmd")
  output_file = root * ".html"

  weave(input_file; doctype = "md2html", out_path = joinpath(Paths.GRAPH, "html"))

  final_file = joinpath(Paths.GRAPH, "html", output_file)

  if Sys.isapple()
    run(`open $final_file`)
  elseif Sys.islinux()
    run(`xdg-open $final_file`)
  else
    @warn "Automatic opening not supported on this OS. File saved at $final_file"
  end
end

####################################################################################################

function list_reports()
  dir = joinpath(Paths.SRC, "live")
  files = filter(f -> endswith(f, ".jmd"), readdir(dir))
  names = String[]
  descriptions = String[]

  for f in files
    root = replace(f, ".jmd" => "")
    push!(names, root)

    lines = readlines(joinpath(dir, f))
    desc_line = findfirst(l -> startswith(strip(l), "# Description"), lines)

    if desc_line !== nothing && desc_line < length(lines)
      desc = strip(lines[desc_line+1])
      push!(descriptions, desc)
    else
      push!(descriptions, "")
    end
  end

  df = DataFrame(name = names, description = descriptions)

  pretty_table(df; display_size = (20, 200))
end

####################################################################################################
