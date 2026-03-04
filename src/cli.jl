####################################################################################################
# src/cli.jl
####################################################################################################

module CLI

####################################################################################################

using ArgParse
using Avicenna.Workflow
using ..Process
using ..Workflows

####################################################################################################

function run_sigma(args)
  s = ArgParseSettings()
  @add_arg_table s begin
    "--config", "-c"
    help = "Path to TOML configuration"
    required = true
    "--no-cache"
    help = "Disable caching"
    action = :store_true
    "--plot"
    help = "Generate plots"
    action = :store_true
  end
  parsed = parse_args(args, s)

  config_dict = Core.load_toml_config(parsed["config"])  # you implement
  cache = Workflow.Cache("cache/sigma", !parsed["no-cache"])

  workflow = parsed["plot"] ? Workflows.sigma_workflow_with_plots : Workflows.sigma_workflow
  result = Workflow.run(workflow, config_dict, cache = cache)

  println("Workflow completed. Cache hits: ", join(result.cache_hits, ", "))
  return 0
end

####################################################################################################

function plot_matrix(args)
  s = ArgParseSettings()
  @add_arg_table s begin
    "--cost"
    required = true

    "--ids"
    required = true

    "--meta"
    required = true

    "--out"
    required = true

    "--pad"
    arg_type = Int, default = 25

    "--title"
    default = ""

  end
  parsed = parse_args(args, s)

  cost = readdlm(parsed["cost"], ',', Float64)
  ids_df = Core.readdf(parsed["ids"]; sep = ',')
  ids = Core.dataframe_to_ids(ids_df)
  meta = Core.readdf(parsed["meta"]; sep = ',')

  plt = Core.plot_grouped_costmatrix(
    cost,
    ids,
    meta;
    pad = parsed["pad"],
    title = parsed["title"],
  )
  savefig(plt, parsed["out"])
  return 0
end

####################################################################################################

function main(args)
  if length(args) == 0
    println("Usage: julia cli.jl <command> [options]")
    println("Commands: run, plot")
    return 1
  end
  command = args[1]
  rest = args[2:end]
  if command == "run"
    run_sigma(rest)
  elseif command == "plot"
    plot_matrix(rest)
  else
    println("Unknown command: $command")
    return 1
  end
end

####################################################################################################

end

####################################################################################################
