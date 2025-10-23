using Clustering
using DataFrames
using Distances
using Plots
using Random
using RDatasets
using UMAP

prostate = dataset("gamair", "prostate")
labels_true = prostate.type
X = Matrix{Float64}(prostate[:, Not([:Column1, :type])])

Random.seed!(42)
embedding = umap(X'; n_neighbors = 15, min_dist = 0.1, metric = Euclidean())

k = length(unique(labels_true))
R = kmeans(embedding, k)
labels_pred = R.assignments

p1 = scatter(
  embedding[1, :],
  embedding[2, :];
  group = labels_true,
  legend = :outertopright,
  title = "UMAP Embedding Colored by True Labels",
  xlabel = "UMAP-1",
  ylabel = "UMAP-2",
)

p2 = scatter(
  embedding[1, :],
  embedding[2, :];
  group = labels_pred,
  legend = :outertopright,
  title = "UMAP Embedding Colored by k-means Clusters",
  xlabel = "UMAP-1",
  ylabel = "UMAP-2",
)

plot(p1, p2, layout = (1, 2), size = (1000, 500))
