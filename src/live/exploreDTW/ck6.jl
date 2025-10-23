let

  using Distances
  using DynamicAxisWarping
  using Plots
  using Statistics

  len_sin = 200
  x_sin = range(0, 2π; length = len_sin)
  y_sin = sin.(x_sin)

  scales = 0.5:0.1:1.5
  cos_lengths = Int.(round.(len_sin .* scales))

  raw_costs = Float64[]
  norm_costs = Float64[]

  for len_cos in cos_lengths
    x_cos = range(0, 2π; length = len_cos)
    y_cos = cos.(x_cos)

    # Use dtw, not dtw_cost (supports unequal lengths)
    cost_val, i1, i2 = dtw(y_sin, y_cos, SqEuclidean())
    push!(raw_costs, cost_val)

    # Normalize by average sequence length
    avg_len = mean([len_sin, len_cos])
    push!(norm_costs, cost_val / avg_len)

    # Show heatmaps for smallest and largest
    if len_cos == minimum(cos_lengths) || len_cos == maximum(cos_lengths)
      D = dtw_cost_matrix(y_sin, y_cos, SqEuclidean())
      plt_heat = heatmap(
        D;
        xlabel = "sin(x) index",
        ylabel = "cos(x) index",
        title = "DTW Cost Matrix (cos length=$len_cos)",
        colorbar_title = "Cost",
      )
      display(plt_heat)
    end
  end

  # Scatterplot of raw DTW cost vs cosine sequence length
  plt_raw = scatter(
    cos_lengths,
    raw_costs;
    xlabel = "Cosine sequence length",
    ylabel = "Raw DTW cost",
    title = "Raw DTW Cost vs Cosine Sequence Length",
    markersize = 6,
    label = "Raw DTW cost",
  )
  plot!(plt_raw, cos_lengths, raw_costs; seriestype = :line, lw = 2, label = false)
  display(plt_raw)

  # Scatterplot of normalized DTW cost vs cosine sequence length
  plt_norm = scatter(
    cos_lengths,
    norm_costs;
    xlabel = "Cosine sequence length",
    ylabel = "Normalized DTW cost",
    title = "Normalized DTW Cost vs Cosine Sequence Length",
    markersize = 6,
    label = "Normalized DTW cost",
  )
  plot!(plt_norm, cos_lengths, norm_costs; seriestype = :line, lw = 2, label = false)
  display(plt_norm)

end
