# # interactive
# include("src/live/saturatedSubsampling/config.jl")
# include("src/live/saturatedSubsampling/util.jl")

# Load path definitions
include(joinpath((pwd() == @__DIR__) ? "src" : "../..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

# Load configuration structs
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))

@info "paths loaded"

using Clustering
using Dates
using DataFrames
using Distances
using DynamicAxisWarping
using LinearAlgebra
using Plots
using Statistics
using Random

@info "dependencies loaded"

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = TrajectoryParams(
  metadata = Vars.SMETA_xlsx,
  batches = [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
)

# Load and split
bundles = load_experiments(sigma_params)

# Collect all metadata into one DataFrame
meta = vcat([b.metadata for b in values(bundles)]...)

# Keep only the columns of interest
meta = select(meta, [:Animal_nr, :Sex, :Genotype])

# Drop rows where Animal_nr == 0 (if those are placeholders)
filter!(row -> row.Animal_nr != 0, meta)
rename!(meta, :Animal_nr => :Animal)

subdfs = split_by_animal(bundles)

# Choose variables
vars = names(subdfs[1])[2:end] .|> Symbol
# vars = [:RT_RER]

# Collect subsamples for all variables at once
subsample_results = collect_subsamples(subdfs, vars, sigma_params)

# sigma_test = [7, 8, 21, 22]
# dfs = Dict(k => subdfs[k] for k in sigma_test)
# subsample_results = collect_subsamples(dfs, vars, sigma_params)

@info "bundle loaded"

# Now run DTW + clustering for each variable
for var in vars

  # var = first(vars)
  println("Processing variable: $var")

  order = sortperm(subsample_results[var].ids; by = id -> (id.subject, id.ixs[1]))
  all_ordered_subsamples = subsample_results[var].subsamples[order]
  all_ordered_ids = subsample_results[var].ids[order]

  N = length(all_ordered_subsamples)
  if N == 0
    @warn "No subsamples collected for $var, skipping"
  end

  # Choose ~10% at random
  k = max(1, round(Int, 0.1 * N))
  idxs = rand(1:N, k)

  # Plot them with fixed y-axis
  plt = plot(
    title = "Random 10% subsamples for $var",
    xlabel = "Time index",
    ylabel = "Value",
    legend = false,
  )

  for i in idxs
    plot!(plt, all_ordered_subsamples[i], alpha = 0.6, lw = 1)
  end
  display(plt)

  # ht = histograms_from_dict(group_pairs(all_ordered_ids))
  # display(ht)

  # Compute DTW cost matrix
  cost_matrix = zeros(Float64, N, N)
  for i = 1:N
    for j = i:N
      cost, _, _ = dtw(all_ordered_subsamples[i], all_ordered_subsamples[j], SqEuclidean())
      norm_cost =
        cost / mean([length(all_ordered_subsamples[i]), length(all_ordered_subsamples[j])])
      cost_matrix[i, j] = norm_cost
      cost_matrix[j, i] = norm_cost
    end
  end

  # save calculations
  writedlm(joinpath(Paths.TMP, string(var, "_cost_matrix.csv")), cost_matrix, ',')
  writedf(joinpath(Paths.TMP, string(var, "_order_ids.csv")), ids_to_dataframe(all_ordered_ids))

  # stats = cost_stats(cost_matrix)
  # print(stats)

  # # Plot cost matrix
  # plt = heatmap(
  #   cost_matrix;
  #   title = "Clustered DTW Costs — $var",
  #   xlabel = "",
  #   ylabel = "",
  #   colorbar_title = "Cost",
  #   size = (800, 700),
  #   xticks = false,
  #   yticks = false,
  #   yflip = true,
  # )
  # display(plt)

  plt = plot_grouped_costmatrix(cost_matrix, all_ordered_ids, meta;)
  display(plt)

  # Hierarchical clustering
  tree = hclust(cost_matrix; linkage = :ward)

  # # Plot cost matrix
  # plt = heatmap(
  #   cost_matrix[tree.order, tree.order];
  #   title = "Clustered DTW Costs — $var",
  #   xlabel = "",
  #   ylabel = "",
  #   colorbar_title = "Cost",
  #   size = (800, 700),
  #   xticks = false,
  #   yticks = false,
  #   yflip = true,
  # )
  # display(plt)

  plt = plot_grouped_costmatrix(
    cost_matrix[tree.order, tree.order],
    all_ordered_ids[tree.order],
    meta;
  )
  display(plt)

end

# # edge detection
# edges_h_ord = diff(cost_matrix_ord; dims = 1)
# edges_v_ord = diff(cost_matrix_ord; dims = 2)

# edges_mag = sqrt.(edges_h_ord[:, 1:end-1] .^ 2 .+ edges_v_ord[1:end-1, :] .^ 2)

# plt = heatmap(edges_mag; color = :inferno, title = "Gradient magnitude (edges)")
# display(plt)

# μ, σ = mean(edges_mag), std(edges_mag)
# mask = edges_mag .> μ + 3σ   # keep only extreme edges
# plt = heatmap(mask; color = :gray, title = "Thresholded edges")
# display(plt)

# row_score = sum(edges_mag, dims = 2)[:]   # boundary strength across rows
# col_score = sum(edges_mag, dims = 1)[:]   # boundary strength across columns

# p1 = plot(
#   row_score;
#   title = "Row boundary score",
#   xlabel = "Index",
#   ylabel = "Score",
#   legend = false,
# )

# p2 = plot(
#   col_score;
#   title = "Column boundary score",
#   xlabel = "Index",
#   ylabel = "Score",
#   legend = false,
# )

# plt = plot(p1, p2; layout = (2, 1), size = (800, 600))
# display(plt)

# plt = heatmap(cost_matrix_ord; color = :inferno, title = "Cost matrix with edges")
# contour!(
#   plt,
#   edges_mag;
#   levels = [μ + 4σ],
#   color = :red,
#   alpha = 0.4,
#   linewidth = 0.5,
#   label = "",
# )
# display(plt)

# # block detection
# A = cost_matrix_ord
# n = size(A, 1)
# w = max(2, round(Int, 0.1 * n))   # 10% window width
# boundary_score = zeros(n - 1)

# # Slide window across columns
# for start = 1:(n-w+1)
#   cols = start:(start+w-1)

#   # Row profile for this window
#   row_profile = mean(A[:, cols]; dims = 2) |> vec

#   # Differences between adjacent rows
#   d = diff(row_profile)

#   # Accumulate absolute differences
#   boundary_score .+= abs.(d)
# end

# # Normalize by number of windows
# boundary_score ./= (n - w + 1)

# # Threshold: μ + 2σ
# μ, σ = mean(boundary_score), std(boundary_score)
# selected = findall(boundary_score .> μ + 4σ)

# # Convert to half-integers for overlay
# boundaries = Float64.(selected) .+ 0.5

# plt = plot(boundary_score; color = :black, lw = 2, label = "score")
# hline!(plt, [μ + 4σ]; color = [:blue], linestyle = :dash, label = ["μ+4σ"])
# display(plt)

# plt = heatmap(A; legend = false, title = "Stability-selected boundaries")
# for b in boundaries
#   vline!(plt, [b]; color = :red, alpha = 0.4, lw = 1, label = false)
#   hline!(plt, [b]; color = :red, alpha = 0.4, lw = 1, label = false)
# end
# display(plt)

# # spectral analysis
# σ = std(cost_matrix)  # scale parameter
# W = exp.(-cost_matrix .^ 2 ./ (2σ^2))  # similarity matrix
# d = sum(W, dims = 2)                # degree vector
# D_inv_sqrt = Diagonal(1 ./ sqrt.(vec(d)))
# Lsym = I - D_inv_sqrt * W * D_inv_sqrt
# using Arpack

# k = 5  # number of clusters
# vals, vecs = eigs(Lsym; nev = k, which = :SM)  # smallest eigenvalues

# X = real(vecs)             # ensure real
# R = kmeans(X', k)          # cluster rows of X
# labels = R.assignments     # cluster assignment per signal
# perm = sortperm(labels)
# plt = heatmap(
#   cost_matrix[perm, perm];
#   title = "Spectral clustering block structure",
#   xlabel = "Index",
#   ylabel = "Index",
# )
# display(plt)
