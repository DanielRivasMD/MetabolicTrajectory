
# row-wise differences
row_mean = [mean(mt[i, :]) for i = 1:size(mt, 1)]
diffs = abs.(diff(row_mean))

# threshold automatically using quantiles
thr = quantile(diffs, 0.99)

edges1 = findall(diffs .> thr)



# hierarchical clustering
hc = hclust(mt; linkage = :average)

labels = cutree(hc; k = 10)

edges2 = findall(diff(labels) .!= 0)



# derivative of cumulative similarity
s = sum(mt, dims = 2)[:]   # row sums
ds = abs.(diff(s))

edges3 = findall(ds .> quantile(ds, 0.99))



using Clustering
using Statistics
using Plots

# ---------------------------------------------------------
# 1. Row‑wise mean differences
# ---------------------------------------------------------
row_mean = [mean(mt[i, :]) for i = 1:size(mt, 1)]
diffs = abs.(diff(row_mean))

thr = quantile(diffs, 0.99)
edges1 = findall(diffs .> thr)

# ---------------------------------------------------------
# 2. Hierarchical clustering
# ---------------------------------------------------------
hc = hclust(mt; linkage = :average)
labels = cutree(hc; k = 5)

edges2 = findall(diff(labels) .!= 0)

# ---------------------------------------------------------
# 3. Row‑sum derivative
# ---------------------------------------------------------
s = sum(mt, dims = 2)[:]
ds = abs.(diff(s))

edges3 = findall(ds .> quantile(ds, 0.99))

# ---------------------------------------------------------
# 4. Comparison Plot
# ---------------------------------------------------------
plt = plot(
  title = "Comparison of Edge Detection Methods",
  xlabel = "Index",
  ylabel = "Signal Value",
  legend = :topright,
  size = (1100, 650),
)

# Plot the three underlying signals
# plot!(plt, row_mean, label="Row Mean", lw=2, color=:blue)
# plot!(plt, s ./ maximum(s), label="Row Sum (scaled)", lw=2, color=:green)
plot!(
  plt,
  labels ./ maximum(labels),
  label = "Cluster Labels (scaled)",
  lw = 2,
  color = :purple,
)

# # Overlay edges from each method
# for e in edges1
#     vline!(plt, [e], color=:red, lw=2, label=false)
# end

for e in edges2
  vline!(plt, [e], color = :orange, lw = 2, label = false)
end

# for e in edges3
#     vline!(plt, [e], color=:black, lw=2, label=false)
# end

# ---------------------------------------------------------
# Legend entries for edges (backend‑safe)
# ---------------------------------------------------------
# scatter!(plt, [0], [0], label="Edges1 (Row Mean)", color=:red, markershape=:none)
scatter!(plt, [0], [0], label = "Edges2 (Clustering)", color = :orange, markershape = :none)
# scatter!(plt, [0], [0], label="Edges3 (Row Sum)", color=:black, markershape=:none)

plt


using Plots

# --- Heatmap of the DTW matrix ---
plt = heatmap(
  mt,
  title = "DTW Matrix with Clustering Boundaries",
  xlabel = "Index",
  ylabel = "Index",
  colorbar = true,
  size = (900, 900),
)

# --- Overlay vertical and horizontal lines for edges2 ---
for e in edges2
  # vertical line
  vline!(plt, [e], color = :orange, lw = 2, label = false)
  # horizontal line
  hline!(plt, [e], color = :orange, lw = 2, label = false)
end

# --- Add a legend entry for the clustering boundaries ---
scatter!(plt, [0], [0], label = "Clustering Edges", color = :orange, markershape = :none)

plt
