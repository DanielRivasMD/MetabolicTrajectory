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
  df = readdf(Vars.SIG1R_HT_csv; sep = ',')
  vars = setdiff(names(df), Vars.xvars_csv)

  for c in Vars.xvars_csv
    df[!, c] = DateTime.(df[!, c], dateformat"mm/dd/yyyy HH:MM:SS")
  end

  for c in vars
    df[!, c] = Float64.(df[!, c])
  end

  # Group by subject
  subdfs = split_by_animal(df)

  for var in vars[1:3]
    println("Processing variable: $var")

    # Normalize variable name (strip trailing _N if present)
    basevar = Symbol(replace(var, r"_\d+$" => ""))

    # Collect subsamples and IDs
    all_subsamples = Vector{Vector{Float64}}()
    all_ids = String[]

    # Iterate over dictionary entries
    for (animal_id, subdf) in subdfs
      # Clean the signal: drop missing values
      raw_signal = subdf[!, basevar]
      signal = collect(skipmissing(raw_signal))

      n = length(signal)
      if n == 0
        @warn "Animal $animal_id has no valid data for $basevar, skipping"
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
        cost, _, _ = dtw(all_subsamples[i], all_subsamples[j], SqEuclidean())
        norm_cost = cost / mean([length(all_subsamples[i]), length(all_subsamples[j])])
        cost_matrix[i, j] = norm_cost
        cost_matrix[j, i] = norm_cost
      end
    end

    # Hierarchical clustering
    tree = hclust(cost_matrix; linkage = :ward)
    leaf_order = tree.order
    cost_matrix_ord = cost_matrix[leaf_order, leaf_order]

    # calcualte edges
    edges_h = diff(cost_matrix; dims = 1)  # horizontal edges (row-to-row changes)
    edges_v = diff(cost_matrix; dims = 2)  # vertical edges (col-to-col changes)

    # --- Plot clustered heatmap ---
    plt = heatmap(
      cost_matrix;
      title = "Clustered DTW Costs — $basevar",
      xlabel = "",
      ylabel = "",
      colorbar_title = "Cost",
      size = (800, 700),
      xticks = false,
      yticks = false,
    )

    plt = heatmap(
      cost_matrix_ord;
      title = "Clustered DTW Costs — $basevar",
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
