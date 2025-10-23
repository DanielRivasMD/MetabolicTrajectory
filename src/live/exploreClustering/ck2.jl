using Clustering
using DataFrames
using Distances
using Plots
using Random
using RDatasets
using TSne

prostate = dataset("gamair", "prostate")
labels_true = prostate.type
X = Matrix{Float64}(prostate[:, Not([:Column1, :type])])

Random.seed!(42)
embedding = tsne(X)'

k = length(unique(labels_true))
R = kmeans(embedding, k)
labels_pred = R.assignments

p1 = scatter(
  embedding[1, :],
  embedding[2, :];
  group = labels_true,
  legend = :outertopright,
  title = "t-SNE Embedding Colored by True Labels",
  xlabel = "t-SNE-1",
  ylabel = "t-SNE-2",
)

p2 = scatter(
  embedding[1, :],
  embedding[2, :];
  group = labels_pred,
  legend = :outertopright,
  title = "t-SNE Embedding Colored by k-means Clusters",
  xlabel = "t-SNE-1",
  ylabel = "t-SNE-2",
)

plot(p1, p2, layout = (1, 2), size = (1000, 500))
