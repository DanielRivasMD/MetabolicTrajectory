#!/usr/bin/env julia

using ArgParse
using DelimitedFiles
using DataFrames
using Clustering
using Colors
using Plots

####################################################################################################
# LOAD CSV AS MATRIX + SUBJECT VECTOR
####################################################################################################

"""
    load_subject_matrix(path::String)

Reads a CSV where:
- Row 1 is: subject,1,2,3,...
- Column 1 is subject IDs
- Remaining entries are Float64

Returns:
- subjects::Vector{Int}
- matrix::Matrix{Float64}
"""
function load_subject_matrix(path::String)
  raw = readdlm(path, ',', Any)

  header = raw[1, 2:end]              # subject IDs in header
  subjects = Int.(header)

  data = raw[2:end, 2:end]            # numeric matrix
  mat = Array{Float64}(data)

  return subjects, mat
end

####################################################################################################
# PLOTTING FUNCTION
####################################################################################################

"""
    plot_subject_heatmap(mat, subjects, meta; pad=20, title="")

Creates a heatmap with:
- core matrix
- top annotation bar for sex
- second annotation bar for genotype
"""
function plot_subject_heatmap(
  mat::Matrix{Float64},
  subjects::Vector{Int},
  meta::DataFrame;
  pad::Int = 20,
  title::String = "",
)
  N = length(subjects)
  @assert size(mat, 1) == N && size(mat, 2) == N

  # Lookups
  sex_lookup = Dict(row.Animal => row.Sex for row in eachrow(meta))
  geno_lookup = Dict(row.Animal => row.Genotype for row in eachrow(meta))

  # Colors
  sex_color(sex) =
    sex == "F" ? RGB(1, 0, 0) : sex == "M" ? RGB(0, 0, 1) : RGB(0.5, 0.5, 0.5)
  geno_color(g) =
    g == "S1RKO" ? RGB(0, 0, 0) :
    g == "Hetero" ? RGB(0.5, 0.5, 0.5) : g == "WT" ? RGB(1, 1, 1) : RGB(0.2, 0.2, 0.2)

  # padded matrix
  core = fill(Float64(NaN), N + pad, N + pad)
  core[pad+1:end, 1:N] .= mat

  pad_colors = fill(RGBA(0, 0, 0, 0), N + pad, N + pad)

  # Two annotation bars: sex + genotype
  block = div(pad, 2)
  sec1 = 1:block
  sec2 = block+1:2block

  for i = 1:N
    s = subjects[i]
    sex = sex_lookup[s]
    geno = geno_lookup[s]

    # Sex bar
    for r in sec1
      pad_colors[r, i] = sex_color(sex)
      pad_colors[pad+i, N+r] = sex_color(sex)
    end

    # Genotype bar
    for r in sec2
      pad_colors[r, i] = geno_color(geno)
      pad_colors[pad+i, N+r] = geno_color(geno)
    end
  end

  plt = heatmap(
    core;
    title = title,
    color = :inferno,
    yflip = true,
    size = (900, 800),
    legend = false,
    nan_color = RGBA(0, 0, 0, 0),
    axis = false,
  )

  plot!(plt, pad_colors; seriestype = :heatmap, yflip = true, legend = false, axis = false)

  return plt
end

####################################################################################################
# MAIN
####################################################################################################

function main()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--csv"
    help = "Input CSV file (e.g. RT_WheelMeters.csv)"
    required = true

    "--meta"
    help = "Metadata CSV (Animal,Sex,Genotype)"
    required = true

    "--out"
    help = "Output PNG path"
    required = true

    "--pad"
    help = "Padding for annotation bars"
    arg_type = Int
    default = 20

    "--title"
    help = "Plot title"
    default = ""
  end

  args = parse_args(s)

  csv_path = args["csv"]
  meta_path = args["meta"]
  out_path = args["out"]
  pad = args["pad"]
  title = args["title"] == "" ? basename(csv_path) : args["title"]

  println("Loading matrix: $csv_path")
  subjects, mat = load_subject_matrix(csv_path)

  println("Loading metadata: $meta_path")
  meta = DataFrame(
    readdlm(meta_path, ',', Any; header = true)[1],
    Symbol.(readdlm(meta_path, ',', Any; header = true)[2]),
  )

  println("Clustering…")
  tree = hclust(mat; linkage = :ward)
  order = tree.order

  mat_ord = mat[order, order]
  subjects_ord = subjects[order]

  println("Plotting…")
  plt = plot_subject_heatmap(mat_ord, subjects_ord, meta; pad = pad, title = title)

  println("Saving → $out_path")
  savefig(plt, out_path)
end

main()

