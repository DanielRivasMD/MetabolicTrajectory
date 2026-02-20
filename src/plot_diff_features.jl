#!/usr/bin/env julia

using DelimitedFiles
using DataFrames
using Plots
using FilePathsBase
using Statistics

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

function readdf(path; sep = '\t')
  data, header = readdlm(path, sep, header = true)
  return DataFrame(data, vec(header))
end

# Convert cumulative vector → per-interval increments
function diff_cumulative(v::Vector{T}) where {T<:Real}
  if length(v) == 0
    return v
  end
  return [v[1]; v[2:end] .- v[1:end-1]]
end

# Plotting logic (separate from looping)
function plot_feature_all_subjects(feature_name, subject_data, subject_colors; outpath)
  # subject_data :: Dict{Int,Vector{Float64}}
  # subject_colors :: Dict{Int,Symbol}

  # ---------------------------------------------------------
  # 1. Pad all subjects to the same length using `missing`
  # ---------------------------------------------------------
  maxlen = maximum(length.(values(subject_data)))

  padded = Dict(
    subj => vcat(vals, fill(missing, maxlen - length(vals))) for
    (subj, vals) in subject_data
  )

  # y-axis scaling from first subject (ignoring missings)
  first_values = skipmissing(first(values(padded)))
  ymax = maximum(first_values) * 1.5

  plt = plot(
    title = feature_name,
    xlabel = "Index",
    ylabel = "Value",
    xlim = (0, 2500),
    ylim = (0, ymax),
    legend = false,
    size = (1000, 600),
  )

  # alpha based on number of subjects
  α = 7 / length(padded)

  # ---------------------------------------------------------
  # 2. Plot each subject
  # ---------------------------------------------------------
  for (subj, vals) in padded
    x = 1:length(vals)
    plot!(plt, x, vals; color = subject_colors[subj], alpha = α, linewidth = 2)
  end

  # ---------------------------------------------------------
  # 3. Compute mean lines (skip missing)
  # ---------------------------------------------------------
  female_vals = [vals for (subj, vals) in padded if subject_colors[subj] == :red]
  male_vals = [vals for (subj, vals) in padded if subject_colors[subj] == :blue]

  if !isempty(female_vals)
    mat = reduce(hcat, female_vals)
    female_mean = [mean(skipmissing(mat[i, :])) for i = 1:maxlen]
    plot!(plt, 1:maxlen, female_mean; color = :red, linewidth = 3, alpha = 1.0)
  end

  if !isempty(male_vals)
    mat = reduce(hcat, male_vals)
    male_mean = [mean(skipmissing(mat[i, :])) for i = 1:maxlen]
    plot!(plt, 1:maxlen, male_mean; color = :blue, linewidth = 3, alpha = 1.0)
  end

  # ---------------------------------------------------------
  # 4. Vertical lines at multiples of 320
  # ---------------------------------------------------------
  for xline = 320:320:2500
    vline!(plt, [xline], color = :black, alpha = 0.3, linewidth = 2)
  end

  mkpath(dirname(outpath))
  savefig(plt, outpath)
end



# -------------------------------------------------------------------
# Main logic
# -------------------------------------------------------------------

function main()
  sigma_dir = "sigma"
  outdir = "graph/png"
  mkpath(outdir)

  # Load metadata
  meta = readdf("sigma/meta.csv"; sep = ',')
  subject_to_color = Dict{Int,Symbol}()

  for row in eachrow(meta)
    subject_to_color[row.Animal] = row.Sex == "F" ? :red : :blue
  end

  # List subject files
  subject_files = filter(f -> occursin(r"exp_\d+\.csv", f), readdir(sigma_dir))

  # Load all subjects into memory
  all_subjects = Dict{Int,DataFrame}()

  for file in subject_files
    subject_id = parse(Int, splitext(split(file, "_")[end])[1])
    df = readdf(joinpath(sigma_dir, file); sep = ',')
    all_subjects[subject_id] = df
  end

  # Extract feature names (columns 2:end) from the first subject
  first_df = first(values(all_subjects))
  feature_names = names(first_df)[2:end]

  # Loop over features
  for feat in feature_names
    println("Plotting feature: $feat")

    # Collect values per subject
    subject_data = Dict{Int,Vector{Float64}}()

    for (subj, df) in all_subjects
      raw_any = df[!, feat]

      # Convert Any → Float64 (handles strings, missing, etc.)
      raw = Float64.(coalesce.(raw_any, 0))

      # cumulative → increments
      vals = diff_cumulative(raw)

      subject_data[subj] = vals
    end

    outpath = joinpath(outdir, "$(feat).png")

    plot_feature_all_subjects(feat, subject_data, subject_to_color; outpath = outpath)
  end

  println("All feature plots saved to graph/png/")
end

main()
