let

  using Distances
  using DynamicAxisWarping
  using Plots
  using Random

  Random.seed!(2025)

  x = range(0, 4π; length = 400)
  signal1 = sin.(x) .+ 0.1 .* randn(length(x))
  signal2 = cos.(x) .+ 0.1 .* randn(length(x))

  radii = [2, 5, 10, 15, 25, 50]

  plots = Plots.Plot[]   # <-- fixed

  for r in radii
    D_r = dtw_cost_matrix(
      signal1,
      signal2,
      SqEuclidean(),
      radiuslimits(r, length(signal1), length(signal2))...,
    )

    plt = heatmap(
      D_r;
      xlabel = "sin(x) index",
      ylabel = "cos(x) index",
      title = "DTW Cost Matrix (radius=$r)",
      colorbar_title = "Cost",
    )
    push!(plots, plt)
  end

  display(plot(plots...; layout = (2, 3), size = (1000, 600)))

end
