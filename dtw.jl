using DataFrames
using DelimitedFiles
using DynamicAxisWarping
using Distances
using Statistics

# ------------------------------------------------------------
# Read CSV using DelimitedFiles
# ------------------------------------------------------------
function read_csv_df(path::String)
  raw = readdlm(path, ',', Any; header = false)
  header = Symbol.(raw[1, :])
  data = raw[2:end, :]
  return DataFrame(data, header)
end

# ------------------------------------------------------------
# Read all exp_*.csv files into Dict{Int,DataFrame}
# ------------------------------------------------------------
function readdf_dict(dir::String)
  dict = Dict{Int,DataFrame}()

  for file in readdir(dir)
    m = match(r"^exp_(\d+)\.csv$", file)
    m === nothing && continue

    key = parse(Int, m.captures[1])
    df = read_csv_df(joinpath(dir, file))
    dict[key] = df
  end

  return dict
end

# ------------------------------------------------------------
# Extract measurement variable names (all columns except Date_Time)
# ------------------------------------------------------------
function measurement_vars(dfs::Dict{Int,DataFrame})
  first_df = first(values(dfs))
  return Symbol.(names(first_df)[2:end])
end

# ------------------------------------------------------------
# Extract a signal vector for a given subject + variable
# ------------------------------------------------------------
extract_signal(df::DataFrame, var::Symbol) = collect(skipmissing(df[!, var]))

# ------------------------------------------------------------
# Build a flat list of all signals across subjects & measurements
# ------------------------------------------------------------
function collect_all_signals(dfs)
  vars = measurement_vars(dfs)

  signals = Vector{Vector{Float64}}()
  labels = Vector{String}()

  for (subj, df) in dfs
    for var in vars
      sig = extract_signal(df, var)
      isempty(sig) && continue

      push!(signals, sig)
      push!(labels, string("subj", subj, "_", var))
    end
  end

  return signals, labels
end

# ------------------------------------------------------------
# Compute full DTW matrix across all signals
# ------------------------------------------------------------
function dtw_full_matrix(signals)
  N = length(signals)
  cost = zeros(Float64, N, N)
  for i = 1:N
    for j = i:N
      s1 = signals[i]
      s2 = signals[j]
      d, _, _ = dtw(s1, s2, SqEuclidean())
      nd = d / mean([length(s1), length(s2)])
      cost[i, j] = nd
      cost[j, i] = nd
    end
  end

  return cost
end

# ------------------------------------------------------------
# Save DTW matrix with labels using DelimitedFiles
# ------------------------------------------------------------
function save_labeled_matrix(path::String, labels, cost)
  N = length(labels)
  out = Array{Any}(undef, N + 1, N + 1)

  out[1, 1] = "signal"

  for j = 1:N
    out[1, j+1] = labels[j]
  end

  for i = 1:N
    out[i+1, 1] = labels[i]
    for j = 1:N
      out[i+1, j+1] = cost[i, j]
    end
  end

  writedlm(path, out, ',')
end

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------
function compute_cross_measurement_dtw(dir::String)
  println("Loading dataframes…")
  dfs = readdf_dict(dir)

  println("Collecting all signals…")
  signals, labels = collect_all_signals(dfs)

  println("Computing full DTW matrix…")
  cost = dtw_full_matrix(signals)

  outdir = joinpath(dir, "full_cross")
  isdir(outdir) || mkdir(outdir)

  outfile = joinpath(outdir, "cross_measurement_dtw.csv")
  println("Saving: $outfile")

  save_labeled_matrix(outfile, labels, cost)

  println("Done.")
end

# ------------------------------------------------------------
# Run it
# ------------------------------------------------------------
sigma_dir = "sigma"
compute_cross_measurement_dtw(sigma_dir)
