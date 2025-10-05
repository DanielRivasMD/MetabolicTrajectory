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

  # convert Any-typed numeric columns to Float64 vectors
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

"plot VO2 and VCO2 over time"
function plot_gas_exchange(df::DataFrame)
  plt = plot(df.DateTime, df.VO2_M, label = "VO₂ (ml/min)", lw = 2)
  plot!(plt, df.DateTime, df.VCO2_M, label = "VCO₂ (ml/min)", lw = 2)
  xlabel!("Time")
  ylabel!("Gas Exchange")
  title!("Gas Exchange Time Series")
  return plt
end

###################################################################################################

"plot locomotor activity (PedMeters_M vs WheelMeters_M)"
function plot_activity(df::DataFrame)
  plt = plot(df.DateTime, df.PedMeters_M, label = "Ped Meters", lw = 2)
  plot!(plt, df.DateTime, df.WheelMeters_M, label = "Wheel Meters", lw = 2)
  xlabel!("Time")
  ylabel!("Distance (m)")
  title!("Locomotor Activity")
  return plt
end

###################################################################################################

"plot environmental variables (Temp, RH, Light)"
function plot_environment(df::DataFrame)
  plt = plot(df.DateTime, df.EnviroTemp_M, label = "Temperature (°C)", lw = 2)
  plot!(plt, df.DateTime, df.EnviroRH_M, label = "Humidity (%)", lw = 2)
  plot!(plt, df.DateTime, df.EnviroLightlux_M, label = "Light (lux)", lw = 2)
  xlabel!("Time")
  ylabel!("Environment")
  title!("Environmental Conditions")
  return plt
end

###################################################################################################

# Example usage (comment out in production)
df = load_timeseries("../data/example.xlsx")
display(plot_gas_exchange(df))
display(plot_activity(df))
display(plot_environment(df))

###################################################################################################
