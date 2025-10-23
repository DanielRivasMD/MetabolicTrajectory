###################################################################################################
# src/DTW.jl
#
# Dynamic Time Warping for metabolic time series data
###################################################################################################

# load config
begin
  include(joinpath("src", "config", "paths.jl"))
  using .Paths
  Paths.ensure_dirs()

  include(joinpath(Paths.CONFIG, "vars.jl"))
  include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
  include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))
  include(joinpath(Paths.UTIL, "DTW.jl"))
end;

###################################################################################################

# load packages
begin
  using DataFrames
  using Dates
  using Plots
  using UnicodePlots
  using DynamicAxisWarping
end;

###################################################################################################

# load dataframe
df = load_timeseries(Vars.SIG1R_HT_file)

# group by subject
gdf = groupby(df, :Animal)

# variables
vars = setdiff(names(df), Vars.exclude_vars)

####################################################################################################


for v in vars
  animals, M = dtw_distance_matrix_by_variable(df, Symbol(v))
  println("\nDTW distances for variable $v")
  display(M)
end

animals, M = dtw_distance_matrix_by_variable(df, :VO2_M)
heatmap(
  1:length(animals),
  1:length(animals),
  M;
  title = "DTW distances for VO2_M",
  xlabel = "Animal index",
  ylabel = "Animal index",
)

####################################################################################################

# DTW analysis per variable
for v in vars
  @info "DTW distance matrix for $(v)"
  animals, M = dtw_distance_matrix_by_variable(df, Symbol(v); window = 30) # window ~ 30 steps band
  plot_dtw_matrix_unicode(animals, M, Symbol(v))
end

####################################################################################################
