###################################################################################################
# src/timeSeries.jl
#
# Visualization of metabolic time series data
###################################################################################################

# load packages
begin
  using DataFrames
  using Dates
  using Plots
  include("../util/ioDataFrame.jl")   # brings in readdf, readxlsx, writedf
end

###################################################################################################

"load and prepare the TimeSeries dataframe"
function load_timeseries(path::AbstractString)
  df = readxlsx(path, sheetname = "TimeSeries")

  # normalize any 1-column matrices into vectors
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
        # skip if conversion fails (e.g. non-numeric column like :Animal)
      end
    end
  end

  # ensure DateTime column is parsed
  if :DateTime in names(df)
    try
      df.DateTime = DateTime.(df.DateTime, dateformat"yyyy-mm-dd HH:MM:SS")
    catch
      df.DateTime = DateTime.(string.(df.DateTime))
    end
  end

  return df
end

###################################################################################################

"plot VO2 and VCO2 over time, grouped by Animal"
function plot_gas_exchange(df::DataFrame)
  plt = plot(legend = false)
  for subdf in groupby(df, :Animal)
    if any(!ismissing, subdf.VO2_M)
      plot!(plt, subdf.DateTime, subdf.VO2_M, lw = 1.5)
    end
    if any(!ismissing, subdf.VCO2_M)
      plot!(plt, subdf.DateTime, subdf.VCO2_M, lw = 1.5, ls = :dash)
    end
  end
  xlabel!("Time")
  ylabel!("Gas Exchange")
  title!("Gas Exchange Time Series by Animal")
  return plt
end

###################################################################################################

"plot locomotor activity (PedMeters_M vs WheelMeters_M), grouped by Animal"
function plot_activity(df::DataFrame)
  plt = plot(legend = false)
  for subdf in groupby(df, :Animal)
    if any(!ismissing, subdf.PedMeters_M)
      plot!(plt, subdf.DateTime, subdf.PedMeters_M, lw = 1.5)
    end
    if any(!ismissing, subdf.WheelMeters_M)
      plot!(plt, subdf.DateTime, subdf.WheelMeters_M, lw = 1.5, ls = :dash)
    end
  end
  xlabel!("Time")
  ylabel!("Distance (m)")
  title!("Locomotor Activity by Animal")
  return plt
end

###################################################################################################

"plot environmental variables (Temp, RH, Light), grouped by Animal"
function plot_environment(df::DataFrame)
  plt = plot(legend = false)
  for subdf in groupby(df, :Animal)
    if any(!ismissing, subdf.EnviroTemp_M)
      plot!(plt, subdf.DateTime, subdf.EnviroTemp_M, lw = 1.5)
    end
    if any(!ismissing, subdf.EnviroRH_M)
      plot!(plt, subdf.DateTime, subdf.EnviroRH_M, lw = 1.5, ls = :dash)
    end
    if any(!ismissing, subdf.EnviroLightlux_M)
      plot!(plt, subdf.DateTime, subdf.EnviroLightlux_M, lw = 1.5, ls = :dot)
    end
  end
  xlabel!("Time")
  ylabel!("Environment")
  title!("Environmental Conditions by Animal")
  return plt
end

###################################################################################################

# Example usage (comment out in production)
df = load_timeseries("../data/example.xlsx")

# ensure graph directory exists
isdir("../graph") || mkpath("../graph")

# generate plots
plt1 = plot_gas_exchange(df)
plt2 = plot_activity(df)
plt3 = plot_environment(df)

# display interactively
display(plt1)
display(plt2)
display(plt3)

# save to files
savefig(plt1, "../graph/gas_exchange.png")
savefig(plt2, "../graph/activity.png")
savefig(plt3, "../graph/environment.png")

###################################################################################################
