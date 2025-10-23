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

  # Load data
  df = readdf(Vars.SIG1R_HT_csv)
  for c in names(df)[2:end]
    df[!, c] = Float64.(df[!, c])
  end


  # Group by subject
  gdf = groupby(df, :Animal)

  # Variables
  vars = setdiff(names(df), Vars.xvars_csv)

  for var in vars[1:3]
    println("Processing variable: $var")

    # Collect subsamples and IDs
    all_subsamples = Vector{Vector{Float64}}()
    all_ids = String[]

    for subdf in gdf
      animal_id = Int(unique(subdf.Animal)[1])

      # Clean the signal: drop missing values
      raw_signal = subdf[:, var]
      signal = collect(skipmissing(raw_signal))

      n = length(signal)
      if n == 0
        @warn "Animal $animal_id has no valid data for $var, skipping"
        continue
      end

      subsample_len = round(Int, 0.05 * n)
      n_samples = 100

      for i = 1:n_samples
        len_var = round(Int, subsample_len * (1 .+ (rand() * 0.02 .- 0.01)))
        start_idx = rand(1:(n-len_var))
        subsample = signal[start_idx:start_idx+len_var-1]

        push!(all_subsamples, subsample)
        push!(all_ids, string(animal_id, "_", start_idx))
      end
    end

    N = length(all_subsamples)

    # Compute DTW cost matrix
    cost_matrix = zeros(Float64, N, N)
    for i = 1:N
      for j = i:N
        cost, i1, i2 = dtw(all_subsamples[i], all_subsamples[j], SqEuclidean())
        norm_cost = cost / mean([length(all_subsamples[i]), length(all_subsamples[j])])
        cost_matrix[i, j] = norm_cost
        cost_matrix[j, i] = norm_cost
      end
    end

    # Hierarchical clustering
    tree = hclust(cost_matrix; linkage = :ward)
    leaf_order = tree.order
    cost_matrix_ord = cost_matrix[leaf_order, leaf_order]

    # --- Plot clustered heatmap ---
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
