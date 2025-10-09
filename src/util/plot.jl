####################################################################################################

# TODO: make multiple dispatch for dataframe & groupeddataframe
function scatter_by_animal(df::DataFrame, var::Symbol)
  gdf = groupby(df, :Animal)
  plt = nothing
  t0 = minimum(df.DateTime)

  for subdf in gdf
    # drop rows with missing in this variable
    sdf = dropmissing(subdf, [var])
    isempty(sdf) && continue

    # convert DateTime to hours since start
    x = (Dates.value.(sdf.DateTime) .- Dates.value(t0)) ./ 3.6e12

    # ensure numeric, skip if conversion fails
    y = try
      Float64.(sdf[!, var])
    catch
      continue
    end

    label = string(first(sdf.Animal))

    if plt === nothing
      plt = scatterplot(
        x,
        y;
        title = string(var),
        xlabel = "Time (h)",
        ylabel = string(var),
        name = label,
        height = displaysize()[1],
        width = displaysize()[2],
      )
    else
      scatterplot!(plt, x, y; name = label)
    end
  end
  return plt
end

####################################################################################################

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


