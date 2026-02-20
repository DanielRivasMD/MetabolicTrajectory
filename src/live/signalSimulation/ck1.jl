####################################################################################################

# Load path definitions
include(joinpath((pwd() == @__DIR__) ? "src" : "../..", "config", "paths.jl"))
using .Paths
Paths.ensure_dirs()

# Load configuration structs
include(joinpath(Paths.CONFIG, "vars.jl"))
include(joinpath(Paths.UTIL, "params.jl"))
include(joinpath(Paths.UTIL, "wrangle.jl"))
include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))

using Clustering
using Dates
using DataFrames
using Distances
using DynamicAxisWarping
using LinearAlgebra
using Plots
using Statistics
using Random

# Sigma experiment: one metadata XLSX, three CSV batches
sigma_params = TrajectoryParams(
  metadata = Vars.SMETA_xlsx,
  batches = [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
  nsamples = 100
)

# Load and split
bundles = load_experiments(sigma_params)

# collect metadata
# Collect all metadata into one DataFrame
meta = vcat([b.metadata for b in values(bundles)]...)

# Keep only the columns of interest
meta = select(meta, [:Animal_nr, :Sex, :Genotype])

# Drop rows where Animal_nr == 0 (if those are placeholders)
filter!(row -> row.Animal_nr != 0, meta)
rename!(meta, :Animal_nr => :Animal)
meta.Group = string.(meta.Sex, "_", meta.Genotype)
animal_to_group = Dict(row.Animal => row.Group for row in eachrow(meta))

subdfs = split_by_animal(bundles)

####################################################################################################

struct SignalParams
  duration::Float64                  # total seconds
  step::Float64                      # time step (1/fs)
  freq_range::Tuple{Float64,Float64} # (min Hz, max Hz)
  amp_range::Tuple{Float64,Float64}  # (min, max amplitude)
  phase_range::Tuple{Float64,Float64}# (min, max phase in radians)
end

function sine_wave(
  amplitude::Float64,
  frequency::Float64,
  phase::Float64,
  duration::Float64,
  step::Float64,
)
  t = 0:step:duration
  ω = 2π * frequency
  y = amplitude .* sin.(ω .* t .+ phase)
  return t, y
end

function build_signal(sp::SignalParams, seed::Int, num_components::Int)
  Random.seed!(seed)
  t = 0:sp.step:sp.duration
  y_total = zeros(Float64, length(t))

  components = []
  for _ = 1:num_components
    f = rand() * (sp.freq_range[2] - sp.freq_range[1]) + sp.freq_range[1]
    a = rand() * (sp.amp_range[2] - sp.amp_range[1]) + sp.amp_range[1]
    p = rand() * (sp.phase_range[2] - sp.phase_range[1]) + sp.phase_range[1]
    push!(components, (frequency = f, amplitude = a, phase = p))
    _, y = sine_wave(a, f, p, sp.duration, sp.step)
    y_total .+= y
  end

  return t, y_total, components
end

function fake_times(n; step_seconds = 270)
  start = DateTime(2024, 1, 1, 0, 0, 0)
  return [start + Dates.Second(step_seconds * (i - 1)) for i = 1:n]
end

####################################################################################################

# --- Experiment setup ---
# Half-day = 12 hours = 43200 seconds
# Each index = 4.5 minutes = 270 seconds
duration_half = 12 * 60 * 60   # 12 hours in seconds
step = 270.0          # 4.5 minutes per index

# Define ranges
sp = SignalParams(duration_half, step, (0.0001, 0.001), (1.0, 5.0), (0.0, 2π))

# First half-day (seed 111)
t1, y1, comps1 = build_signal(sp, 111, 10)
t2, y2, comps2 = build_signal(sp, 222, 10)
day_signal1 = vcat(y1, y2)
week_signal1 = repeat(day_signal1, 7)

# Second week with different seeds
t3, y3, comps3 = build_signal(sp, 333, 10)
t4, y4, comps4 = build_signal(sp, 444, 10)
day_signal2 = vcat(y3, y4)
week_signal2 = repeat(day_signal2, 7)

t5, y5, comps5 = build_signal(sp, 555, 1)
t6, y6, comps6 = build_signal(sp, 666, 1)
day_signal3 = vcat(y5, y6)
week_signal3 = repeat(day_signal3, 7)

# Plot one week
plt = plot(
  1:(Int(length(week_signal1) / 7)),
  week_signal1[1:(Int(length(week_signal1) / 7))],
  title = "Week Signal 1",
  legend = false,
)
display(plt)
plt = plot(
  1:(Int(length(week_signal2) / 7)),
  week_signal2[1:(Int(length(week_signal1) / 7))],
  title = "Week Signal 2",
  legend = false,
)
display(plt)
plt = plot(
  1:(Int(length(week_signal3) / 7)),
  week_signal3[1:(Int(length(week_signal1) / 7))],
  title = "Week Signal 3",
  legend = false,
)
display(plt)

plt = plot(1:length(week_signal1), week_signal1, title = "Week Signal 1", legend = false)
display(plt)
plt = plot(1:length(week_signal2), week_signal2, title = "Week Signal 2", legend = false)
display(plt)
plt = plot(1:length(week_signal3), week_signal3, title = "Week Signal 3", legend = false)
display(plt)

####################################################################################################

s1 = collect_subsamples(week_signal1, fake_times(length(week_signal1)), sigma_params)
s2 = collect_subsamples(week_signal2, fake_times(length(week_signal2)), sigma_params)
s3 = collect_subsamples(week_signal3, fake_times(length(week_signal3)), sigma_params)

merged = merge_subsamplecontainers([s1, s2, s3], [1, 2, 3])

order = sortperm(merged.ids; by = id -> (id.subject, id.ixs[1]))
ordered_subsamples = merged.subsamples[order]
ordered_ids = merged.ids[order]

N = length(ordered_subsamples)
if N == 0
  @warn "No subsamples collected for $var, skipping"
end

# Compute DTW cost matrix
cost_matrix = zeros(Float64, N, N)

for i = 1:N
  for j = i:N
    cost, _, _ = dtw(ordered_subsamples[i], ordered_subsamples[j], SqEuclidean())
    norm_cost = cost / mean([length(ordered_subsamples[i]), length(ordered_subsamples[j])])
    cost_matrix[i, j] = norm_cost
    cost_matrix[j, i] = norm_cost
  end
end

plt = plot_grouped_costmatrix(cost_matrix, ordered_ids, meta)
display(plt)

# Hierarchical clustering
tree = hclust(cost_matrix; linkage = :ward)

plt = plot_grouped_costmatrix(cost_matrix[tree.order, tree.order], ordered_ids[tree.order], meta)
display(plt)

####################################################################################################

day_signal1 = cumsum(day_signal1)
day_signal2 = cumsum(day_signal2)
day_signal3 = cumsum(day_signal3)

week_signal1 = repeat(day_signal1, 7)
week_signal2 = repeat(day_signal2, 7)
week_signal3 = repeat(day_signal3, 7)

# Plot one week
plt = plot(
  1:(Int(length(week_signal1) / 7)),
  week_signal1[1:(Int(length(week_signal1) / 7))],
  title = "Week Signal 1",
  legend = false,
)
display(plt)
plt = plot(
  1:(Int(length(week_signal2) / 7)),
  week_signal2[1:(Int(length(week_signal1) / 7))],
  title = "Week Signal 2",
  legend = false,
)
display(plt)
plt = plot(
  1:(Int(length(week_signal3) / 7)),
  week_signal3[1:(Int(length(week_signal1) / 7))],
  title = "Week Signal 3",
  legend = false,
)
display(plt)

plt = plot(1:length(week_signal1), week_signal1, title = "Week Signal 1", legend = false)
display(plt)
plt = plot(1:length(week_signal2), week_signal2, title = "Week Signal 2", legend = false)
display(plt)
plt = plot(1:length(week_signal3), week_signal3, title = "Week Signal 3", legend = false)
display(plt)

s1 = collect_subsamples(week_signal1, fake_times(length(week_signal1)), sigma_params)
s2 = collect_subsamples(week_signal2, fake_times(length(week_signal2)), sigma_params)
s3 = collect_subsamples(week_signal3, fake_times(length(week_signal3)), sigma_params)

merged = merge_subsamplecontainers([s1, s2, s3], [1, 2, 3])

order = sortperm(merged.ids; by = id -> (id.subject, id.ixs[1]))
ordered_subsamples = merged.subsamples[order]
ordered_ids = merged.ids[order]

N = length(ordered_subsamples)
if N == 0
  @warn "No subsamples collected for $var, skipping"
end

# Compute DTW cost matrix
cost_matrix = zeros(Float64, N, N)

for i = 1:N
  for j = i:N
    cost, _, _ = dtw(ordered_subsamples[i], ordered_subsamples[j], SqEuclidean())
    norm_cost = cost / mean([length(ordered_subsamples[i]), length(ordered_subsamples[j])])
    cost_matrix[i, j] = norm_cost
    cost_matrix[j, i] = norm_cost
  end
end

plt = plot_grouped_costmatrix(cost_matrix, ordered_ids, meta)
display(plt)

# Hierarchical clustering
tree = hclust(cost_matrix; linkage = :ward)

plt = plot_grouped_costmatrix(cost_matrix[tree.order, tree.order], ordered_ids[tree.order], meta)
display(plt)

####################################################################################################
