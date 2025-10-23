###################################################################################################
# src/timeSeries.jl
#
# Fast Fourier Transform for metabolic time series data
###################################################################################################

# load config
begin
  include(joinpath("src", "config", "paths.jl"))
  using .Paths
  Paths.ensure_dirs()

  include(joinpath(Paths.CONFIG, "vars.jl"))
  include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
  include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))
  include(joinpath(Paths.UTIL, "FFT.jl"))
end;

###################################################################################################

# load packages
begin
  using DataFrames
  using Dates
  using Plots
  using UnicodePlots
  using FFTW
end;

###################################################################################################

# load dataframe
df = load_timeseries(Vars.SIG1R_HT_file)

# group by subject
gdf = groupby(df, :Animal)

# variables
vars = setdiff(names(df), Vars.exclude_vars)

####################################################################################################

# TODO: hardcoded
fs = 0.00056

for v in vars
  spectra = fft_by_animal(df, Symbol(v); fs = fs)
  @info "FFT for $v"
  for (animal, (freqs, mags)) in spectra
    println("Animal $animal, variable $v, first 5 freqs: ", freqs[1:5])
    println("Magnitudes: ", mags[1:5])
  end
end

for v in vars
  spectra = fft_by_animal(df, Symbol(v); fs = fs)
  @info "FFT for $v"

  for (animal, (freqs, mags)) in spectra
    # only keep the first half (positive frequencies)
    N = length(freqs)
    half = 1:div(N, 2)

    f = freqs[half]
    m = mags[half] ./ N   # normalize magnitude

    println("Animal $animal, variable $v")
    display(
      lineplot(
        f,
        m;
        xlabel = "Frequency (Hz)",
        ylabel = "Magnitude",
        title = "FFT of $v (Animal $animal)",
      ),
    )
  end
end

####################################################################################################
