let

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

  using Clustering
  using Dates
  using DataFrames
  using Distances
  using DynamicAxisWarping
  using LinearAlgebra
  using Plots
  using Statistics
  using Random
  using UMAP

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
  meta.Group = string.(meta.Sex, "_", meta.Genotype)
  animal_to_group = Dict(row.Animal => row.Group for row in eachrow(meta))

  subdfs = split_by_animal(bundles)

  # # DOC: hardcoded sigma test
  # sigma_test = [7, 8, 21, 22]

  # Choose variables
  # DOC: hardcoded variables
  vars = names(subdfs[1])[2:end] .|> Symbol
  # vars = [:RT_RER]

  # Collect subsamples for all variables at once
  subsample_results = collect_subsamples(subdfs, vars, sigma_params)
  # dfs = Dict(k => subdfs[k] for k in sigma_test)
  # subsample_results = collect_subsamples(dfs, vars, sigma_params)

  # Now run DTW + clustering for each variable
  for var in vars

    # var = first(vars)
    println("Processing variable: $var")

    order = sortperm(subsample_results[var].ids; by = split_id)
    all_ordered_subsamples = subsample_results[var].subsamples[order]
    all_ordered_ids = subsample_results[var].ids[order]
    prefixes = parse.(Int, first.(split.(all_ordered_ids, "_")))
    # sufixes = parse.(Int, last.(split.(all_ordered_ids, "_")))
    groups = [animal_to_group[p] for p in prefixes]

    N = length(all_ordered_subsamples)
    if N == 0
      @warn "No subsamples collected for $var, skipping"
    end

    # # Choose ~1% at random
    # k = max(1, round(Int, 0.01 * N))
    # idxs = rand(1:N, k)

    # # Plot them with fixed y-axis
    # plt = plot(
    #   title = "Random 1% subsamples for $var",
    #   xlabel = "Time index",
    #   ylabel = "Value",
    #   legend = false,
    # )

    # for i in idxs
    #   plot!(plt, all_ordered_subsamples[i], alpha = 0.6, lw = 1)
    # end
    # display(plt)

    # ht = histograms_from_dict(group_pairs(all_ordered_ids))
    # display(ht)

    # Compute DTW cost matrix
    cost_matrix = zeros(Float64, N, N)
    for i = 1:N
      for j = i:N
        cost, _, _ =
          dtw(all_ordered_subsamples[i], all_ordered_subsamples[j], SqEuclidean())
        norm_cost =
          cost /
          mean([length(all_ordered_subsamples[i]), length(all_ordered_subsamples[j])])
        cost_matrix[i, j] = norm_cost
        cost_matrix[j, i] = norm_cost
      end
    end

    using DelimitedFiles
    using FilePathsBase: mkpath

    # ensure output directory exists
    outdir = joinpath(pwd(), "csv")
    mkpath(outdir)

    # build filename using the variable name
    outfile = joinpath(outdir, string(var, "_cost_matrix.csv"))

    # write the matrix
    writedlm(outfile, cost_matrix, ',')
    println("Saved cost matrix for $var → $outfile")

    stats = cost_stats(cost_matrix)
    print(stats)

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
    # )
    # display(plt)

    # plt = plot_grouped_costmatrix(cost_matrix, groups, repeat(collect(1:100), 4))
    # display(plt)

    # # Hierarchical clustering
    # tree = hclust(cost_matrix; linkage = :ward)
    # leaf_order = tree.order
    # cost_matrix_ord = cost_matrix[leaf_order, leaf_order]
    # groups_ord = groups[leaf_order]

    # # Plot clustered heatmap
    # plt = heatmap(
    #   cost_matrix_ord;
    #   title = "Clustered DTW Costs — $var",
    #   xlabel = "",
    #   ylabel = "",
    #   colorbar_title = "Cost",
    #   size = (800, 700),
    #   xticks = false,
    #   yticks = false,
    # )
    # display(plt)

    # plt = plot_grouped_costmatrix(cost_matrix_ord, groups_ord, leaf_order)
    # display(plt)

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

  # # umap
  # # Inputs
  # C = cost_matrix            # n × n, symmetric, nonnegative
  # n = size(C, 1)

  # # UMAP on precomputed distances
  # Random.seed!(42)
  # embedding = umap(
  #   C;
  #   n_neighbors = 15,
  #   min_dist = 0.1,
  #   metric = :precomputed,   # key: treat C as distances
  # )  # returns 2 × n (or n × 2 depending on your UMAP; adjust indexing accordingly)

  # # Choose number of clusters (tune as needed)
  # k = 5
  # R = kmeans(embedding, k)     # k-means on rows (points)
  # labels_pred = R.assignments

  # # Build plots
  # p1 = scatter(
  #   embedding[1, :],
  #   embedding[2, :];
  #   group = labels_pred,
  #   legend = :outertopright,
  #   title = "UMAP embedding colored by k-means clusters",
  #   xlabel = "UMAP-1",
  #   ylabel = "UMAP-2",
  # )

  # p2 = scatter(
  #   embedding[1, :],
  #   embedding[2, :];
  #   group = 1:n,                 # proxy: index-based coloring
  #   legend = :none,
  #   title = "UMAP embedding colored by signal index",
  #   xlabel = "UMAP-1",
  #   ylabel = "UMAP-2",
  # )

  # plt = plot(p1, p2; layout = (1, 2), size = (1000, 500))
  # display(plt)

  # # PCA clustering from DTW costs
  # pca_res = pca_clustering_from_cost(cost_matrix_ord; n_components = 10, k = 5)
  # perm_pca = sortperm(pca_res.assignments)
  # plt = heatmap(
  #   cost_matrix_ord[perm_pca, perm_pca];
  #   title = "PCA clustering block structure",
  #   xlabel = "",
  #   ylabel = "",
  #   legend = false,
  # )
  # display(plt)

  # # t-SNE clustering from DTW costs
  # tsne_res = tsne_clustering_from_cost(cost_matrix_ord; perplexity = 30, k = 5)

  # # Reorder indices by cluster assignment
  # perm_tsne = sortperm(tsne_res.assignments)

  # # Visualize t-SNE embedding
  # plt_embed = scatter(
  #   tsne_res.embedding[:, 1],
  #   tsne_res.embedding[:, 2];
  #   group = tsne_res.assignments,
  #   legend = :outertopright,
  #   title = "t-SNE embedding colored by k-means clusters",
  #   xlabel = "t-SNE-1",
  #   ylabel = "t-SNE-2",
  # )
  # display(plt_embed)

  # # Reordered cost matrix by t-SNE clusters
  # plt_tsne = heatmap(
  #   cost_matrix_ord[perm_tsne, perm_tsne];
  #   title = "t-SNE clustering block structure",
  #   xlabel = "",
  #   ylabel = "",
  #   legend = false,
  # )
  # display(plt_tsne)

end
