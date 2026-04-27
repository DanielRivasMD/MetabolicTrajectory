####################################################################################################

module HCCore

####################################################################################################

using DataFrames
using DelimitedFiles
using Clustering
using Plots
using ..USCore: readdf

####################################################################################################

export cluster_analysis

####################################################################################################

"""
    cluster_analysis(matrix_path, meta_path, output_dir; k=2, plot_title="") -> Dict

Loads the DTW matrix and metadata, performs hierarchical clustering (Ward linkage),
saves a dendrogram, a clustered heatmap, and a CSV with cluster assignments.
"""
function cluster_analysis(
  matrix_path::String,
  meta_path::String,
  output_dir::String;
  k::Int = 2,
  plot_title::String = "",
)
  mkpath(output_dir)

  # 1. Load data
  matrix, subjects = load_dtw_matrix(matrix_path)
  meta = readdf(meta_path; sep = ',')

  # 2. Hierarchical clustering
  hc = hclust(matrix; linkage = :ward)

  # 3. Cut tree into clusters
  labels = cutree(hc; k = k)
  # labels is a vector of length N, indices follow original order

  # 4. Save cluster assignments
  assignments = DataFrame(subject = subjects, cluster = labels)
  assign_path = joinpath(output_dir, "cluster_assignments.csv")
  writedf(assign_path, assignments; sep = ',')

  # 5. Reorder matrix according to clustering and save heatmap
  order = hc.order
  matrix_clustered = matrix[order, order]
  subjects_ordered = subjects[order]

  # Minimal heatmap (no annotation bars, just matrix + meta annotations?) 
  # We'll generate a simple heatmap and then overlay sex/genotype colours manually.
  plt_heat =
    plot_clustered_heatmap(matrix_clustered, subjects_ordered, meta; title = plot_title)
  savefig(plt_heat, joinpath(output_dir, "clustered_heatmap.png"))

  # 6. Dendrogram
  plt_dendro = plot(hc; legend = false, title = plot_title)
  savefig(plt_dendro, joinpath(output_dir, "dendrogram.png"))

  return Dict(
    "output_dir" => output_dir,
    "assignments" => assign_path,
    "dendrogram" => joinpath(output_dir, "dendrogram.png"),
    "heatmap" => joinpath(output_dir, "clustered_heatmap.png"),
    "n_subjects" => length(subjects),
    "k" => k,
  )
end

####################################################################################################

function load_dtw_matrix(path::String)
  raw = readdlm(path, ',', Any)
  header = raw[1, :]
  subjects = Int.(header[2:end])   # first cell is "subject"
  data = Float64.(raw[2:end, 2:end])
  return data, subjects
end

####################################################################################################

function plot_clustered_heatmap(
  matrix::Matrix{Float64},
  subjects::Vector{Int},
  meta::DataFrame;
  title::String = "",
)
  # Build subject lookup
  sex_dict = Dict(row.Animal => row.Sex for row in eachrow(meta))
  geno_dict = Dict(row.Animal => row.Genotype for row in eachrow(meta))

  # Map subjects to colors
  sex_colors(s) = s == "F" ? :red : s == "M" ? :blue : :gray
  geno_colors(g) =
    g == "WT" ? :green : g == "S1RKO" ? :black : g == "Hetero" ? :purple : :gray

  N = length(subjects)
  pad = 20
  core = fill(NaN, N + pad, N + pad)
  core[pad+1:end, 1:N] .= matrix
  core[pad+1:end, N+1:end] .= NaN   # right pad

  # Annotation rows/cols
  ann_rows = fill(RGBA(0, 0, 0, 0), N + pad, N + pad)
  block = div(pad, 2)
  sec1 = 1:block
  sec2 = block+1:2block
  for i = 1:N
    subj = subjects[i]
    sex = sex_dict[subj]
    geno = geno_dict[subj]
    for r in sec1
      ann_rows[r, i] = sex_colors(sex)
      ann_rows[pad+i, N+r] = sex_colors(sex)
    end
    for r in sec2
      ann_rows[r, i] = geno_colors(geno)
      ann_rows[pad+i, N+r] = geno_colors(geno)
    end
  end

  plt = heatmap(
    core;
    color = :inferno,
    title = title,
    yflip = true,
    size = (800, 700),
    legend = false,
    nan_color = RGBA(0, 0, 0, 0),
    axis = false,
  )
  plot!(plt, ann_rows; seriestype = :heatmap, yflip = true, legend = false, axis = false)
  return plt
end

####################################################################################################

end

####################################################################################################
