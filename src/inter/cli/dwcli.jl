####################################################################################################

module DWCLI

####################################################################################################

using ArgParse
using Avicenna.Flow: Cache, launch
using ..DWFlow: flow

####################################################################################################

export run

####################################################################################################

"""
    run(args::Vector{String}) -> Int

CLI entry point for the Dynamic Time Warping workflow.
Usage: avicenna dw --input-dir <path> [--output-dir <path>] [--variables var1,var2] [--no-cache]
"""
function run(args::Vector{String})
  s = ArgParseSettings()
  @add_arg_table! s begin
    "--input-dir", "-i"
    help = "Directory containing per‑animal CSVs and meta.csv"
    required = true
    "--output-dir", "-o"
    help = "Output directory for DTW matrices"
    default = "dtw_full"
    "--variables"
    help = "Comma‑separated list of variables to process (default: all)"
    default = nothing
    "--no-cache"
    help = "Disable caching"
    action = :store_true
  end
  parsed = parse_args(args, s)

  input_dir = parsed["input-dir"]
  if !isdir(input_dir)
    @error "Input directory not found: $input_dir"
    return 1
  end

  # Split variables if provided
  var_list = if parsed["variables"] !== nothing
    split(parsed["variables"], ",") .|> strip
  else
    nothing
  end

  flow_config =
    Dict{String,Any}("input_dir" => input_dir, "output_dir" => parsed["output-dir"])
  if var_list !== nothing
    flow_config["variables"] = var_list
  end

  cache = Cache("cache/dynamic_time_warping", !parsed["no-cache"])
  result = launch(flow, flow_config; cache = cache)

  println("DTW analysis completed.")
  summary = result.stage_outputs["02_compute_dtw"]
  println("Processed variables: ", join(summary["variables_processed"], ", "))
  println("Output directory: ", summary["output_dir"])
  return 0
end

####################################################################################################

end

####################################################################################################
