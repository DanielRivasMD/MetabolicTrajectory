####################################################################################################

# src/document.jl
####################################################################################################

module Document

####################################################################################################

using Weave
using Avicenna.Workflow
using ..Process

####################################################################################################

function report(result::Workflow.WorkflowResult, outpath::String; format = "html")
  # Extract information from result.stage_outputs
  dtw_results = result.stage_outputs["dtw"]
  # Build Markdown with summaries and perhaps include plots
  md = """
  # Sigma Analysis Report

  ## Configuration
  ```
  $(result.origin["config"])
  ```

  ## Results per variable
  """
  for (var, data) in dtw_results
    md *= "\n### $var\n"
    md *= "- Matrix size: $(size(data.matrix,1))×$(size(data.matrix,2))\n"
    # Optionally embed a plot (Weave can handle file paths)
  end

  tmpfile = tempname() * ".jmd"
  write(tmpfile, md)

  if format == "html"
    weave(tmpfile; doctype = "md2html", out_path = outpath)
  elseif format == "pdf"
    weave(tmpfile; doctype = "md2pdf", out_path = outpath)
  end
end

####################################################################################################

end

####################################################################################################
