###################################################################################################

using DataFrames
using DelimitedFiles
using XLSX

###################################################################################################

"read dataframe"
function readdf(path; sep = '\t')
  data, header = readdlm(path, sep, header = true)
  return DataFrame(data, vec(header))
end

###################################################################################################

"read dataframe"
function readxlsx(path; sheetname::AbstractString = "TimeSeries")
  xf = XLSX.readxlsx(path)
  sheet = xf[sheetname]

  # extract header row and data block
  header = vec(sheet[1, :]) |> collect
  data = sheet[2:end, :] |> collect

  df = DataFrame(data, Symbol.(header))  # use symbols for column names
  return df
end

###################################################################################################

"write dataframe"
function writedf(path, df::DataFrame; sep = ',')
  header = permutedims(names(df))  # 1×N matrix of strings
  data = Matrix(df)
  writedlm(path, vcat(header, data), sep)
end

###################################################################################################
