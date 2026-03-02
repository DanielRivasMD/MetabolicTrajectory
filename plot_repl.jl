using DelimitedFiles
using DataFrames
using Clustering
using Plots
using StatsPlots
using Random

include("src/util/ioDataFrame.jl")
include("src/util/params.jl")
include("src/util/wrangle.jl")

# Set random seed for reproducibility
Random.seed!(42)

# File paths
base_path = "dtw/raw/n100_len80.000_var0.050_limits641_1920"
cost_path = joinpath(base_path, "RT_WheelMeters_cost_matrix.csv")
ids_path = joinpath(base_path, "RT_WheelMeters_order_ids.csv")
meta_path = "sigma/meta.csv"

println("Loading cost matrix: $cost_path")
cost_matrix = Array{Float64}(readdlm(cost_path, ',', Float64))

println("Loading order IDs: $ids_path")
all_ordered_ids = dataframe_to_ids(readdf(ids_path; sep = ','))

println("Loading metadata: $meta_path")
meta = readdf(meta_path; sep = ',')

println("Performing hierarchical clustering…")
tree = hclust(cost_matrix; linkage = :ward)

println("Clustering order indices: ", tree.order)

# Reorder the matrix according to clustering
cost_matrix_ordered = cost_matrix[tree.order, tree.order]
ids_ordered = all_ordered_ids[tree.order]

# Set number of clusters
n_clusters = 4
cluster_labels = cutree(tree; k = n_clusters)

# Print cluster information
println("\n=== CLUSTER INFORMATION ===")
for c = 1:n_clusters
  cluster_indices = findall(cluster_labels .== c)
  cluster_samples = tree.order[cluster_indices]
  println("Cluster $c: $(length(cluster_indices)) samples")
  println("  Sample positions: $(cluster_samples[1:min(5, end)])...")
  println()
end

# Create color palette for clusters
cluster_colors = [:red, :blue, :green, :purple, :orange, :brown, :pink, :cyan][1:n_clusters]

# -------------------------------------------------------------------
# 1. DENDOGRAM (colored leaves)
# -------------------------------------------------------------------
println("\nCreating dendrogram...")
leaf_colors = cluster_colors[cluster_labels[tree.order]]

dendrogram_plt = plot(
  tree,
  color = leaf_colors,
  linecolor = :gray,
  linewidth = 1.5,
  title = "Dendrogram with $n_clusters Clusters",
  size = (1200, 400),
  legend = false,
  margin = 5Plots.mm,
)

# Add cluster legend
for c = 1:n_clusters
  n_samples = count(cluster_labels .== c)
  scatter!(
    [0],
    [0],
    color = cluster_colors[c],
    markershape = :circle,
    markersize = 5,
    label = "Cluster $c (n=$n_samples)",
  )
end

display(dendrogram_plt)

# -------------------------------------------------------------------
# 2. HEATMAP (with yflip=true to match dendrogram orientation)
# -------------------------------------------------------------------
println("\nCreating heatmap...")

# The heatmap with yflip=true means row 1 is at TOP
# The dendrogram's tree.order puts first item at TOP of dendrogram
# So they correspond directly when both use yflip=true

heatmap_plt = heatmap(
  cost_matrix_ordered,
  title = "Clustered Cost Matrix (yflip=true)",
  color = :inferno,
  xlabel = "Column Index (dendrogram order)",
  ylabel = "Row Index (dendrogram order)",
  size = (800, 700),
  yflip = true,  # This makes row 1 at TOP
  aspect_ratio = :equal,
)

display(heatmap_plt)

# -------------------------------------------------------------------
# 3. SAMPLE 5 RANDOM SUBSAMPLES FROM EACH CLUSTER
# -------------------------------------------------------------------
println("\nSampling 5 random subsamples from each cluster...")

# Store all plots for the 5x1 layout
cluster_plots = []

# Function to extract signal from original files
function extract_signal(subject::Int, start_idx::Int, end_idx::Int)
  # Path to original subject file
  subject_path = joinpath("sigma", "exp_$(subject).csv")

  if !isfile(subject_path)
    @warn "Subject file not found: $subject_path"
    return nothing
  end

  # Read the subject's data
  df = readdf(subject_path; sep = ',')

  # Extract WheelMeters column (assuming it's the 2nd column or named appropriately)
  # Adjust column name if different
  if "WheelMeters" in names(df)
    signal_raw = df[!, "WheelMeters"]
  elseif size(df, 2) >= 2
    signal_raw = df[!, 2]  # Use second column
  else
    @warn "Cannot find WheelMeters column in $subject_path"
    return nothing
  end

  # Convert to Float64 and handle missing
  signal = Float64.(coalesce.(signal_raw, 0))

  # Differentiate if needed (since we're in diff/ directory)
  # signal = diff_cumulative(signal)  # Uncomment if signals need differentiation

  # Extract the subsample
  if start_idx <= length(signal) && end_idx <= length(signal)
    return signal[start_idx:end_idx]
  else
    @warn "Indices out of bounds for subject $subject: $start_idx-$end_idx (signal length=$(length(signal)))"
    return nothing
  end
end

# Process each cluster
for c = 1:n_clusters
  # Find all indices in this cluster (in original order, not tree.order)
  cluster_indices_original = findall(cluster_labels .== c)

  # Convert to indices in tree.order (for ids_ordered)
  cluster_indices_ordered =
    [findfirst(==(idx), tree.order) for idx in cluster_indices_original]

  # Randomly sample up to 5 from this cluster
  n_to_sample = min(5, length(cluster_indices_ordered))
  sampled_positions = rand(cluster_indices_ordered, n_to_sample)

  println("\nCluster $c: sampling $n_to_sample subsamples")

  # Create a subplot for this cluster
  subplot = plot(
    title = "C-$c (S-$(join(unique([ids_ordered[i].subject for i in sampled_positions]), ',')))",
    xlabel = "Time Index",
    ylabel = "Wheel Meters",
    legend = false,
    size = (800, 200),
  )

  # Extract and plot each sampled subsample
  for (i, pos) in enumerate(sampled_positions)
    id = ids_ordered[pos]
    subject = id.subject
    start_idx, end_idx = id.ixs

    println("  Sampling: Subject $subject, indices $start_idx-$end_idx")

    signal = extract_signal(subject, start_idx, end_idx)

    if !isnothing(signal)
      # Create x-axis (relative time)
      x_vals = 1:length(signal)

      # Use color intensity based on sample (different shades of cluster color)
      alpha_val = 0.3 + 0.7 * (i / n_to_sample)
      plot!(
        subplot,
        x_vals,
        signal,
        color = cluster_colors[c],
        linewidth = 2,
        alpha = alpha_val,
        label = false,
      )
    end
  end

  push!(cluster_plots, subplot)
end

# Combine all cluster plots into a 5x1 layout
cluster_plots_combined =
  plot(cluster_plots..., layout = (5, 1), size = (1000, 1200), margin = 5Plots.mm)

display(cluster_plots_combined)

# -------------------------------------------------------------------
# 4. COMBINE ALL VISUALIZATIONS
# -------------------------------------------------------------------
println("\nCombining all plots...")

# Top row: dendrogram
# Middle: heatmap
# Bottom: cluster samples (5 rows)

combined_plt = plot(
  dendrogram_plt,
  heatmap_plt,
  cluster_plots_combined,
  layout = grid(3, 1, heights = [0.25, 0.35, 0.4]),
  size = (1200, 1800),
  margin = 5Plots.mm,
  title = "Complete Analysis - $n_clusters Clusters",
)

display(combined_plt)

# -------------------------------------------------------------------
# 5. SAVE ALL FIGURES
# -------------------------------------------------------------------
println("\nSaving figures...")

# Create output directory
out_dir = "analysis_output"
mkpath(out_dir)

savefig(dendrogram_plt, joinpath(out_dir, "01_dendrogram.png"))
savefig(heatmap_plt, joinpath(out_dir, "02_heatmap.png"))
savefig(cluster_plots_combined, joinpath(out_dir, "03_cluster_samples.png"))
savefig(combined_plt, joinpath(out_dir, "04_complete_analysis.png"))

println("All figures saved to: $out_dir/")

# -------------------------------------------------------------------
# 6. CORRESPONDENCE BETWEEN HEATMAP AND DENDOGRAM
# -------------------------------------------------------------------
println("\n" * "="^60)
println("CORRESPONDENCE BETWEEN VISUALIZATIONS:")
println("="^60)
println("""
- The dendrogram shows hierarchical clustering with colored leaves
- The heatmap uses yflip=true (row 1 at TOP, row N at BOTTOM)
- The dendrogram's tree.order lists items from TOP to BOTTOM
- Therefore, the TOP of dendrogram corresponds to TOP row of heatmap
- The cluster samples show actual signals from random subsamples
""")

# Print the mapping for verification
println("\nFirst 10 mappings (dendrogram order → original index):")
for i = 1:min(10, length(tree.order))
  orig_idx = tree.order[i]
  cluster = cluster_labels[orig_idx]
  subject = ids_ordered[i].subject
  println("  Dendro pos $i → original idx $orig_idx → Cluster $cluster → Subject $subject")
end

println("\nDone!")
