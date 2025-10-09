###################################################################################################
# src/timeSeries.jl
#
# Visualization of metabolic time series data
###################################################################################################

# load config
begin
  include(joinpath("src", "config", "paths.jl"))
  using .Paths
  Paths.ensure_dirs()

  include(joinpath(Paths.CONFIG, "vars.jl"))
  include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
  include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))
end;

###################################################################################################

# load packages
begin
  using DataFrames
  using Dates
  using Plots
  using UnicodePlots
end;

###################################################################################################

df = load_timeseries(Vars.SIG1R_HT_file)

gdf = groupby(df, :Animal)


using Dates, UnicodePlots, DataFrames

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
      )
    else
      scatterplot!(plt, x, y; name = label)
    end
  end
  return plt
end

vars = setdiff(names(df), ["DateTime", "DurationMin", "Animal"])

for v in vars
  @info v
  plt = scatter_by_animal(df, Symbol(v))
  display(plt)
end


###################################################################################################

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
