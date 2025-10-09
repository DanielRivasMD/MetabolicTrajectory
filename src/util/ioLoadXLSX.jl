####################################################################################################

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
end;

####################################################################################################
