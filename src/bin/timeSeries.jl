###################################################################################################
# src/timeSeries.jl
#
# Visualization of metabolic time series data
###################################################################################################

begin
  # Load path definitions
  include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
  using .Paths
  Paths.ensure_dirs()

  # Load configuration structs
  include(joinpath(Paths.CONFIG, "vars.jl"))
  include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
  include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))
  include(joinpath(Paths.UTIL, "plot.jl"))
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

# load dataframe
df = load_timeseries(Vars.SIG1R_HT_file)

# group by subject
gdf = groupby(df, :Animal)

# variables
vars = setdiff(names(df), Vars.exclude_vars)

# explor variables
for v in vars
  @info v
  plt = scatter_by_animal(df, Symbol(v))
  display(plt)
end

###################################################################################################

# collect stats per subject
for subdf in gdf
  println("Animal = ", first(subdf.Animal))
  display(describe(subdf[:, vars]))
end

###################################################################################################
