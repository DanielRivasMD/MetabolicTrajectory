####################################################################################################

using ArgParse
using DelimitedFiles
using DataFrames
using Clustering
using Plots

####################################################################################################

include("src/util/params.jl")
include("src/util/wrangle.jl")
include("src/util/ioDataFrame.jl")

####################################################################################################

function main()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--cost"
    help = "Path to cost matrix CSV"
    required = true

    "--ids"
    help = "Path to order IDs CSV"
    required = true

    "--meta"
    help = "Path to metadata CSV"
    required = true

    "--out"
    help = "Output PNG file"
    required = true

    "--pad"
    help = "Padding for plot"
    arg_type = Int
    default = 200

    "--title"
    help = "Plot title"
    default = ""
  end

  args = parse_args(s)

  cost_path = args["cost"]
  ids_path = args["ids"]
  meta_path = args["meta"]
  out_path = args["out"]
  pad = args["pad"]
  title = args["title"] == "" ? basename(cost_path) : args["title"]

  println("Loading cost matrix: $cost_path")
  cost_matrix = Array{Float64}(readdlm(cost_path, ',', Float64))

  println("Loading order IDs: $ids_path")
  all_ordered_ids = dataframe_to_ids(readdf(ids_path; sep = ','))

  println("Loading metadata: $meta_path")
  meta = readdf(meta_path; sep = ',')

  println("Performing hierarchical clustering…")
  tree = hclust(cost_matrix; linkage = :ward)

  println("Plotting…")
  plt = plot_grouped_costmatrix(
    cost_matrix[tree.order, tree.order],
    all_ordered_ids[tree.order],
    meta;
    pad = pad,
    title = title,
  )

  println("Saving PNG → $out_path")
  savefig(plt, out_path)
end

####################################################################################################

main()

####################################################################################################
