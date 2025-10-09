###################################################################################################
# src/selfCompareMatrix.jl
#
# Build a self-comparison matrix for VO2_M values and render as heatmap
###################################################################################################

begin
  using DataFrames
  using Plots
  include("../util/ioDataFrame.jl")
end

###################################################################################################

"load and prepare the TimeSeries dataframe"
function load_timeseries(path::AbstractString)
  df = readxlsx(path, sheetname = "TimeSeries")

  # normalize 1-column matrices
  for c in names(df)
    if df[!, c] isa AbstractMatrix
      df[!, c] = vec(df[!, c])
    end
  end

  # replace "." with missing and convert to Float64 where possible
  for c in names(df)
    if eltype(df[!, c]) == Any
      df[!, c] = replace(df[!, c], "." => missing)
      try
        df[!, c] = Float64.(df[!, c])
      catch
        # skip non-numeric columns
      end
    end
  end

  return df
end

###################################################################################################

"build a self-comparison matrix for a numeric vector"
function self_compare_matrix(x::Vector{<:Real})
  n = length(x)
  mat = Array{Float64}(undef, n, n)
  for i = 1:n
    for j = 1:n
      # direct comparison: difference, product, or equality
      # here we use product as the comparison
      mat[i, j] = x[i] * x[j]
    end
  end
  return mat
end

###################################################################################################

"plot the self-comparison matrix as a heatmap"
function plot_self_compare(x::Vector{<:Real}; title_str = "Self-comparison Matrix")
  mat = self_compare_matrix(x)
  plt = heatmap(
    1:length(x),
    1:length(x),
    mat;
    xlabel = "Index",
    ylabel = "Index",
    c = :viridis,
    colorbar_title = "Value",
    title = title_str,
    legend = false,
    aspect_ratio = 1,
  )
  return plt
end

###################################################################################################

# Example usage
df = load_timeseries("../data/example.xlsx")

# take VO2_M column, drop missings
x = collect(skipmissing(df.VO2_M))

isdir("../graph") || mkpath("../graph")

plt = plot_self_compare(x; title_str = "VO₂ Self-comparison Matrix")
display(plt)
savefig(plt, "../graph/vo2_self_compare.png")

###################################################################################################
