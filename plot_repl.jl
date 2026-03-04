using DelimitedFiles
using DataFrames
using Clustering
using Plots
using StatsPlots
using Random
using Statistics

include("src/util/ioDataFrame.jl")
include("src/util/params.jl")
include("src/util/wrangle.jl")

# Convert cumulative vector → per-interval increments
function diff_cumulative(v::Vector{T}) where {T<:Real}
  if length(v) == 0
    return v
  end
  return [v[1]; v[2:end] .- v[1:end-1]]
end

# Set random seed for reproducibility
Random.seed!(42)

# File paths
base_path = "dtw/diff/n100_len320.000_var0.050_limits641_1920"
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

# Create color palette for 4 clusters
cluster_colors = [:red, :blue, :green, :purple]

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

heatmap_plt = heatmap(
  cost_matrix_ordered,
  title = "Clustered Cost Matrix (yflip=true)",
  color = :inferno,
  xlabel = "Column Index (dendrogram order)",
  ylabel = "Row Index (dendrogram order)",
  size = (800, 700),
  yflip = true,
  aspect_ratio = :equal,
)

display(heatmap_plt)

# -------------------------------------------------------------------
# 3. MEAN ± STD FOR EACH CLUSTER
# -------------------------------------------------------------------
println("\nComputing mean and standard deviation for each cluster...")

# Store all plots for the 4x1 layout
cluster_plots = []

# Function to extract and differentiate signal from original files
function extract_diff_signal(subject::Int, start_idx::Int, end_idx::Int)
  # Path to original subject file
  subject_path = joinpath("sigma", "exp_$(subject).csv")

  if !isfile(subject_path)
    @warn "Subject file not found: $subject_path"
    return nothing
  end

  # Read the subject's data
  df = readdf(subject_path; sep = ',')

  # Extract WheelMeters column
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

  # Apply diff_cumulative to get per-interval increments
  signal_diff = diff_cumulative(signal)

  # Extract the subsample
  if start_idx <= length(signal_diff) && end_idx <= length(signal_diff)
    return signal_diff[start_idx:end_idx]
  else
    @warn "Indices out of bounds for subject $subject: $start_idx-$end_idx (signal length=$(length(signal_diff)))"
    return nothing
  end
end

# Process each cluster
for c = 1:n_clusters
  println("\nProcessing Cluster $c...")

  # Find all indices in this cluster (in tree.order order)
  cluster_indices_original = findall(cluster_labels .== c)

  # Get the corresponding positions in ids_ordered
  cluster_positions = [findfirst(==(idx), tree.order) for idx in cluster_indices_original]

  println("  Found $(length(cluster_positions)) subsamples in cluster $c")

  # Collect all signals for this cluster
  cluster_signals = []
  valid_lengths = Int[]

  for pos in cluster_positions
    id = ids_ordered[pos]
    subject = id.subject
    start_idx, end_idx = id.ixs

    signal_diff = extract_diff_signal(subject, start_idx, end_idx)

    if !isnothing(signal_diff)
      push!(cluster_signals, signal_diff)
      push!(valid_lengths, length(signal_diff))
    end
  end

  if isempty(cluster_signals)
    println("  WARNING: No valid signals found for cluster $c")
    continue
  end

  # Determine common length for alignment (use minimum length)
  common_length = minimum(valid_lengths)
  println("  Using common length: $common_length (from min of signals)")

  # Trim all signals to common length
  trimmed_signals = [sig[1:common_length] for sig in cluster_signals]

  # Create matrix: rows = time points, columns = samples
  signal_matrix = hcat(trimmed_signals...)

  # Compute mean and std across samples for each time point
  mean_signal = mean(signal_matrix, dims = 2)[:]
  std_signal = std(signal_matrix, dims = 2)[:]

  # Get unique subjects for title
  unique_subjects = unique([ids_ordered[pos].subject for pos in cluster_positions])
  subject_str = join(string.("S-", unique_subjects[1:min(5, end)]), ", ")
  if length(unique_subjects) > 5
    subject_str *= " … ($(length(unique_subjects)) total)"
  end

  # Create x-axis
  x_vals = 1:common_length

  # Create subplot for this cluster
  subplot = plot(
    title = "C-$c: $(length(cluster_positions)) subsamples ($subject_str)",
    xlabel = "Time Index",
    ylabel = "Δ Wheel Meters",
    legend = :topright,
    size = (800, 200),
  )

  # Plot individual signals with low alpha (background)
  for sig in trimmed_signals
    plot!(
      subplot,
      x_vals,
      sig,
      color = cluster_colors[c],
      linewidth = 1,
      alpha = 0.15,
      label = false,
    )
  end

  # Plot mean as bold line
  plot!(
    subplot,
    x_vals,
    mean_signal,
    color = cluster_colors[c],
    linewidth = 4,
    label = "Mean",
  )

  # Plot mean ± std as shaded region
  plot!(
    subplot,
    x_vals,
    mean_signal + std_signal,
    fillrange = mean_signal - std_signal,
    fillalpha = 0.3,
    linealpha = 0.0,
    color = cluster_colors[c],
    label = "Mean ± Std",
  )

  push!(cluster_plots, subplot)
end

# Combine all cluster plots into a 4x1 layout
cluster_plots_combined = plot(
  cluster_plots...,
  layout = (4, 1),
  size = (1000, 1200),  # Taller to accommodate more info
  margin = 5Plots.mm,
)

display(cluster_plots_combined)

# -------------------------------------------------------------------
# 4. COMBINE ALL VISUALIZATIONS
# -------------------------------------------------------------------
println("\nCombining all plots...")

# Top row: dendrogram
# Middle: heatmap
# Bottom: cluster mean plots (4 rows)

combined_plt = plot(
  dendrogram_plt,
  heatmap_plt,
  cluster_plots_combined,
  layout = grid(3, 1, heights = [0.2, 0.3, 0.5]),  # More space for mean plots
  size = (1200, 1800),
  margin = 5Plots.mm,
  title = "Complete Analysis - $n_clusters Clusters (Mean ± Std)",
)

display(combined_plt)

# -------------------------------------------------------------------
# 5. SAVE ALL FIGURES
# -------------------------------------------------------------------
println("\nSaving figures...")

# Create output directory
out_dir = "analysis_output_mean_std"
mkpath(out_dir)

savefig(dendrogram_plt, joinpath(out_dir, "01_dendrogram.png"))
savefig(heatmap_plt, joinpath(out_dir, "02_heatmap.png"))
savefig(cluster_plots_combined, joinpath(out_dir, "03_cluster_mean_std.png"))
savefig(combined_plt, joinpath(out_dir, "04_complete_analysis_mean_std.png"))

println("All figures saved to: $out_dir/")

# -------------------------------------------------------------------
# 6. SUMMARY STATISTICS
# -------------------------------------------------------------------
println("\n" * "="^60)
println("CLUSTER SUMMARY STATISTICS:")
println("="^60)

for c = 1:n_clusters
  cluster_indices = findall(cluster_labels .== c)
  n_samples = length(cluster_indices)

  # Get subjects in this cluster
  subjects_in_cluster =
    unique([ids_ordered[findfirst(==(idx), tree.order)].subject for idx in cluster_indices])

  # Get metadata for these subjects
  sex_counts = Dict{String,Int}()
  geno_counts = Dict{String,Int}()

  for s in subjects_in_cluster
    row = meta[meta.Animal.==s, :]
    if nrow(row) > 0
      sex = row[1, :Sex]
      geno = row[1, :Genotype]
      sex_counts[sex] = get(sex_counts, sex, 0) + 1
      geno_counts[geno] = get(geno_counts, geno, 0) + 1
    end
  end

  println("Cluster $c: $n_samples subsamples from $(length(subjects_in_cluster)) subjects")
  println("  Sex: $(join(["$k=$v" for (k,v) in sex_counts], ", "))")
  println("  Genotype: $(join(["$k=$v" for (k,v) in geno_counts], ", "))")
  println()
end

println("\nDone! Mean ± standard deviation plotted for each cluster.")
