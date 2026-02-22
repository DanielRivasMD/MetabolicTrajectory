
using DataFrames
using CSV
using Distances   # for DTWDistance()
using StatsBase

# ------------------------------------------------------------
# 1. Read all exp_*.csv files into a Dict{Int,DataFrame}
# ------------------------------------------------------------
function readdf_dict(dir::String; sep = ',')
  dict = Dict{Int,DataFrame}()

  for file in readdir(dir)
    m = match(r"^exp_(\d+)\.csv$", file)
    m === nothing && continue

    key = parse(Int, m.captures[1])
    df = CSV.read(joinpath(dir, file), DataFrame; delim = sep)
    dict[key] = df
  end

  return dict
end

# ------------------------------------------------------------
# 2. Extract variable names (all columns except Date_Time)
# ------------------------------------------------------------
function measurement_vars(dfs::Dict{Int,DataFrame})
  first_df = first(values(dfs))
  return Symbol.(names(first_df)[2:end])
end

# ------------------------------------------------------------
# 3. Extract a vector for a given subject + variable
# ------------------------------------------------------------
function extract_signal(df::DataFrame, var::Symbol)
  return collect(skipmissing(df[!, var]))
end

# ------------------------------------------------------------
# 4. Compute DTW cost matrix for one measurement across subjects
# ------------------------------------------------------------
function dtw_cost_matrix(signals::Dict{Int,Vector{Float64}})
  subjects = collect(keys(signals))
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
# 5. Main function: compute DTW for all variables
# ------------------------------------------------------------
function compute_all_dtw(dir::String)
  println("Loading dataframes…")
  dfs = readdf_dict(dir)

  println("Extracting measurement variables…")
  vars = measurement_vars(dfs)

  dtw_results = Dict{Symbol,Dict}()

  for var in vars
    println("Processing variable: $var")

    # Extract signals for all subjects
    signals = Dict{Int,Vector{Float64}}()
    for (id, df) in dfs
      sig = extract_signal(df, var)
      isempty(sig) && continue
      signals[id] = sig
    end

    # Compute DTW cost matrix
    subjects, cost = dtw_cost_matrix(signals)

    dtw_results[var] = Dict(:subjects => subjects, :cost_matrix => cost)
  end

  return dtw_results
end

# ------------------------------------------------------------
# 6. Run it
# ------------------------------------------------------------
sigma_dir = "sigma"
dtw_results = compute_all_dtw(sigma_dir)
println("DTW computation complete.")
