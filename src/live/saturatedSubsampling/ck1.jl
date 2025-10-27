let

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

  # collect metadata
  # Collect all metadata into one DataFrame
  meta = vcat([b.metadata for b in values(bundles)]...)

  # Keep only the columns of interest
  meta = select(meta, [:Animal_nr, :Sex, :Genotype])

  # Drop rows where Animal_nr == 0 (if those are placeholders)
  filter!(row -> row.Animal_nr != 0, meta)
  rename!(meta, :Animal_nr => :Animal)


  subdfs = split_by_animal(bundles)

  # Choose variables
  # DOC: hardcoded variables
  vars = [:VO2_1, :VCO2_1]

  # Collect subsamples for all variables at once
  subsample_results = collect_subsamples(subdfs, vars, sigma_params)

  # Now run DTW + clustering for each variable
  for var in vars
    println("Processing variable: $var")

    all_subsamples = subsample_results[var].subsamples
    all_ids = subsample_results[var].ids

    N = length(all_subsamples)
    if N == 0
      @warn "No subsamples collected for $var, skipping"
      continue
    end

    # Choose ~1% at random
    k = max(1, round(Int, 0.01 * N))
    idxs = rand(1:N, k)

    # Plot them with fixed y-axis
    plt = plot(
      title = "Random 1% subsamples for $var",
      xlabel = "Time index",
      ylabel = "Value",
      legend = false,
      ylims = (0, 10),
    )

    for i in idxs
      plot!(plt, all_subsamples[i], alpha = 0.6, lw = 1)
    end
    display(plt)

    # Compute DTW cost matrix
    cost_matrix = zeros(Float64, N, N)
    for i = 1:N
      for j = i:N
        cost, _, _ = dtw(all_subsamples[i], all_subsamples[j], SqEuclidean())
        norm_cost = cost / mean([length(all_subsamples[i]), length(all_subsamples[j])])
        cost_matrix[i, j] = norm_cost
        cost_matrix[j, i] = norm_cost
      end
    end

    stats = cost_stats(cost_matrix)
    print(stats)

    # Plot cost matrix
    plt = heatmap(
      cost_matrix;
      title = "Clustered DTW Costs — $var",
      xlabel = "",
      ylabel = "",
      colorbar_title = "Cost",
      size = (800, 700),
      xticks = false,
      yticks = false,
    )
    display(plt)

    # Hierarchical clustering
    tree = hclust(cost_matrix; linkage = :ward)
    leaf_order = tree.order
    cost_matrix_ord = cost_matrix[leaf_order, leaf_order]

    # Plot clustered heatmap
    plt = heatmap(
      cost_matrix_ord;
      title = "Clustered DTW Costs — $var",
      xlabel = "",
      ylabel = "",
      colorbar_title = "Cost",
      size = (800, 700),
      xticks = false,
      yticks = false,
    )
    display(plt)
  end

end
