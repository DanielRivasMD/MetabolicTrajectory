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

    "--sex"
    help = "Filter by sex (M or F)"
    default = ""

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
  sex_filter = args["sex"]
  title = args["title"] == "" ? basename(cost_path) : args["title"]

  println("Loading cost matrix: $cost_path")
  cost_matrix = Array{Float64}(readdlm(cost_path, ',', Float64))

  println("Loading order IDs: $ids_path")
  all_ordered_ids = dataframe_to_ids(readdf(ids_path; sep = ','))

  println("Loading metadata: $meta_path")
  meta = readdf(meta_path; sep = ',')

  if sex_filter != ""
    println("Filtering by sex = $sex_filter")

    # subjects with the desired sex
    allowed_subjects = meta.Animal[meta.Sex.==sex_filter]

    # indices of subsamples whose subject is allowed
    keep = findall(id -> id.subject in allowed_subjects, all_ordered_ids)

    # filter cost matrix and ids
    cost_matrix = cost_matrix[keep, keep]
    all_ordered_ids = all_ordered_ids[keep]
  end

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
