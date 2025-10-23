using Weave
using DataFrames
using PrettyTables

# include paths.jl relative to this file
include(joinpath(@__DIR__, "..", "config", "paths.jl"))
using .Paths

"""
    report(root::String)

Given the root name of an analysis (without extension),
this function weaves `Paths.SRC/live/<root>.jmd` into
`Paths.GRAPH/html/<root>.html` and opens it in the default browser.
Supports both Ubuntu/Linux and macOS.
"""
function report(root::String)
  input_file = joinpath(Paths.SRC, "live", root * ".jmd")
  output_file = root * ".html"

  # Weave to HTML
  weave(input_file; doctype = "md2html", out_path = joinpath(Paths.GRAPH, "html"))

  # Construct full output path
  final_file = joinpath(Paths.GRAPH, "html", output_file)

  # Open in default browser depending on OS
  if Sys.isapple()
    run(`open $final_file`)
  elseif Sys.islinux()
    run(`xdg-open $final_file`)
  else
    @warn "Automatic opening not supported on this OS. File saved at $final_file"
  end
end

"""
    list_reports()

Scans `Paths.SRC/live/` for `.jmd` files and returns a DataFrame
with two columns:
- `:name` — the root name of the report
- `:description` — the text following a `# Description` tag in the file,
  or an empty string if not found.

Also prints the DataFrame using PrettyTables with custom display size.
"""
function list_reports()
  dir = joinpath(Paths.SRC, "live")
  files = filter(f -> endswith(f, ".jmd"), readdir(dir))
  names = String[]
  descriptions = String[]

  for f in files
    root = replace(f, ".jmd" => "")
    push!(names, root)

    # Read file and look for "# Description"
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

  # PrettyTables display with larger display size and no cropping
  pretty_table(df; display_size = (20, 200))
end
