####################################################################################################

using DataFrames, Dates

####################################################################################################

"load and prepare the TimeSeries dataframe"
function load_timeseries(path::AbstractString)
  df = readxlsx(path, sheetname = "TimeSeries")

  df.DateTime = DateTime.(df.DateTime, dateformat"yyyy/mm/dd HH:MM:SS")

  for c in names(df)
    col = df[!, c]

    if c == :Animal
      df[!, c] = string.(col)

      # other string‑typed columns: try to coerce to Float64 with missing
    elseif eltype(col) <: AbstractString
      df[!, c] = map(x -> x == "." ? missing : tryparse(Float64, x), col)

      # Any‑typed columns: attempt the same
    elseif eltype(col) == Any
      df[!, c] = map(x -> x == "." ? missing : tryparse(Float64, string(x)), col)
    end
  end

  return df
end;

####################################################################################################
