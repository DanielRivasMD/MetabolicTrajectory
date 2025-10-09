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
end;

###################################################################################################

# Example usage (comment out in production)
df = load_timeseries(Vars.SIG1R_HT_file)

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
