let

# Seed RNG for reproducibility
Random.seed!(2025)

# Generate noisy sinusoidal signals
x = range(0, 4π; length=400)
signal1 = sin.(x) .+ 0.1 .* randn(length(x))
signal2 = cos.(x) .+ 0.1 .* randn(length(x))

# Plot raw signals
plt_signals = plot(x, signal1; label="noisy sin(x)", color=:blue, lw=2)
plot!(plt_signals, x, signal2; label="noisy cos(x)", color=:red, lw=2,
      title="Noisy Sinusoidal Signals")
display(plt_signals)

# Define radii to test
radii = [2, 5, 10, 15, 25, 50]

# Collect results
results = Dict{Int, Float64}()

for r in radii
    cost_r = dtw_cost(signal1, signal2, SqEuclidean(), r)
    results[r] = cost_r
    println("DTW cost with radius=$r: ", cost_r)
end

# For one example radius, also visualize cost matrix + alignment
chosen_radius = 25
D_r = dtw_cost_matrix(signal1, signal2, SqEuclidean())
plt_heat = heatmap(D_r;
    xlabel="sin(x) index",
    ylabel="cos(x) index",
    title="DTW Cost Matrix (radius=$chosen_radius)",
    colorbar_title="Cost")
display(plt_heat)

cost_val, i1, i2 = dtw(signal1, signal2, SqEuclidean())
plt_align = plot(x, signal1; label="noisy sin(x)", color=:blue, lw=2)
plot!(plt_align, x, signal2; label="noisy cos(x)", color=:red, lw=2)
for (idx1, idx2) in zip(i1[1:15:end], i2[1:15:end])
    plot!(plt_align, [x[idx1], x[idx2]], [signal1[idx1], signal2[idx2]];
          color=:gray, alpha=0.4, label=false)
end
display(plt_align)

# Visualize cost vs radius
plt_costs = plot(radii, [results[r] for r in radii];
    seriestype=:scatter,
    markersize=6,
    xlabel="Radius",
    ylabel="DTW Cost",
    title="DTW Cost vs Radius (Noisy Sinusoids)",
    label="DTW cost")
plot!(plt_costs, radii, [results[r] for r in radii]; seriestype=:line, lw=2, label=false)
display(plt_costs)

end
