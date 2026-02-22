using DataFrames
using DelimitedFiles
using Distances
using Statistics
using Printf

# ------------------------------------------------------------
# 1. Read CSV using DelimitedFiles
# ------------------------------------------------------------
function read_csv_df(path::String)
  raw = readdlm(path, ',', Any; header = false)

  # First row is header
  header = Symbol.(raw[1, :])
  data = raw[2:end, :]

  # Convert to DataFrame
  df = DataFrame(data, header)

  return df
end

# ------------------------------------------------------------
# 2. Read all exp_*.csv files into Dict{Int,DataFrame}
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
# 3. Extract measurement variable names (all columns except Date_Time)
# ------------------------------------------------------------
function measurement_vars(dfs::Dict{Int,DataFrame})
  first_df = first(values(dfs))
  return Symbol.(names(first_df)[2:end])
end

# ------------------------------------------------------------
# 4. Extract a signal vector for a given subject + variable
# ------------------------------------------------------------
extract_signal(df::DataFrame, var::Symbol) = collect(skipmissing(df[!, var]))

# ------------------------------------------------------------
# 5. Compute DTW cost matrix for one measurement across subjects
# ------------------------------------------------------------
function dtw_cost_matrix(signals::Dict{Int,Vector{Float64}})
  subjects = sort(collect(keys(signals)))
  N = length(subjects)

  cost = zeros(Float64, N, N)

  for i = 1:N
    for j = i:N
      s1 = signals[subjects[i]]
      s2 = signals[subjects[j]]
      d = evaluate(DTWDistance(), s1, s2)
      cost[i, j] = d
      cost[j, i] = d
    end
  end

  return subjects, cost
end

# ------------------------------------------------------------
# 6. Save DTW matrix with subject IDs using DelimitedFiles
# ------------------------------------------------------------
function save_dtw_matrix(path::String, var::Symbol, subjects, cost)
  N = length(subjects)

  # Build a matrix with headers
  out = Array{Any}(undef, N + 1, N + 1)

  # Top-left corner
  out[1, 1] = "subject"

  # Column headers
  for j = 1:N
    out[1, j+1] = subjects[j]
  end

  # Row headers + cost matrix
  for i = 1:N
    out[i+1, 1] = subjects[i]
    for j = 1:N
      out[i+1, j+1] = cost[i, j]
    end
  end

  # Write CSV
  writedlm(joinpath(path, string(var) * ".csv"), out, ',')
end

# ------------------------------------------------------------
# 7. Main function: compute DTW for all variables and save results
# ------------------------------------------------------------
function compute_all_dtw(dir::String)
  println("Loading dataframes…")
  dfs = readdf_dict(dir)

  println("Extracting measurement variables…")
  vars = measurement_vars(dfs)

  # Create output directory
  outdir = joinpath(dir, "full")
  isdir(outdir) || mkdir(outdir)

  println("Computing DTW for each variable…")

  for var in vars
    println(" → $var")

    # Extract signals for all subjects
    signals = Dict{Int,Vector{Float64}}()
    for (id, df) in dfs
      sig = extract_signal(df, var)
      isempty(sig) && continue
      signals[id] = sig
    end

    # Compute DTW
    subjects, cost = dtw_cost_matrix(signals)

    # Save CSV
    save_dtw_matrix(outdir, var, subjects, cost)
  end

  println("All DTW matrices saved in: $outdir")
end

# ------------------------------------------------------------
# 8. Run it
# ------------------------------------------------------------
sigma_dir = "sigma"
compute_all_dtw(sigma_dir)
