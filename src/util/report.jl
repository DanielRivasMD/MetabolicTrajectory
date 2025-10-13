using Weave

"""
    report(root::String)

Given the root name of an analysis (without extension),
this function weaves `src/live/<root>.jmd` into
`graph/html/<root>.html` and opens it in the default browser.
Supports both Ubuntu/Linux and macOS.
"""
function report(root::String)
    input_file  = joinpath("src/live", root * ".jmd")
    output_file = root * ".html"

    # Weave to HTML
    weave(input_file;
          doctype="md2html",
          out_path="graph/html")

    # Construct full output path
    final_file = joinpath("graph/html", output_file)

    # Open in default browser depending on OS
    if Sys.isapple()
        run(`open $final_file`)
    elseif Sys.islinux()
        run(`xdg-open $final_file`)
    else
        @warn "Automatic opening not supported on this OS. File saved at $final_file"
    end
end

using DataFrames

"""
    list_reports()

Scans `src/live/` for `.jmd` files and returns a DataFrame
with two columns:
- `:name` — the root name of the report
- `:description` — the text following a `# Description` tag in the file,
  or an empty string if not found.
"""
function list_reports()
    dir = "src/live"
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
            # Take the line immediately after "# Description"
            desc = strip(lines[desc_line + 1])
            push!(descriptions, desc)
        else
            push!(descriptions, "")
        end
    end

    return DataFrame(name = names, description = descriptions)
end
