####################################################################################################

using FFTW, Statistics

####################################################################################################

function fft_by_animal(df::DataFrame, var::Symbol; fs::Float64)
  gdf = groupby(df, :Animal)
  results = Dict()

  for subdf in gdf
    y = subdf[!, var]
    y = skipmissing(y) |> collect
    isempty(y) && continue

    N = length(y)
    Y = fft(y)
    freqs = (0:N-1) .* (fs / N)   # frequency bins

    # store magnitude spectrum
    results[first(subdf.Animal)] = (freqs, abs.(Y))
  end
  return results
end

####################################################################################################
