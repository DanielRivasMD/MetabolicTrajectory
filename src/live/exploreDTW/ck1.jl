let

using Distances
using DynamicAxisWarping
using Plots

# Generate signals
x = range(0, 5π; length=1000)
y_sin = sin.(x)
y_cos = cos.(x)

# Compute DTW cost matrix
D = dtw_cost_matrix(y_sin, y_cos, SqEuclidean())

# Capture heatmap under a variable
plt_heat = heatmap(D;
    xlabel="sin(x) index",
    ylabel="cos(x) index",
    title="DTW Cost Matrix Heatmap",
    colorbar_title="Cost"
)

# Explicitly display it
display(plt_heat)

# Alignment path overlay on original signals
cost, i1, i2 = dtw(y_sin, y_cos, SqEuclidean())
plt_align = plot(x, y_sin; label="sin(x)", color=:blue, lw=2)
plot!(plt_align, x, y_cos; label="cos(x)", color=:red, lw=2)

# Draw alignment lines for a subset
for (idx_sin, idx_cos) in zip(i1[1:20:end], i2[1:20:end])
    plot!(plt_align, [x[idx_sin], x[idx_cos]], [y_sin[idx_sin], y_cos[idx_cos]];
          color=:gray, alpha=0.4, label=false)
end

display(plt_align)

println("Euclidean distance: ", evaluate(Euclidean(), y_sin, y_cos))
println("DTW cost: ", cost)

end
