###################################################################################################

using DataFrames
using DelimitedFiles
using XLSX

###################################################################################################

function readdf(path; sep = '\t')
  data, header = readdlm(path, sep, header = true)
  return DataFrame(data, vec(header))
end

###################################################################################################

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

function writedf(path, df::DataFrame; sep = ',')
  header = permutedims(names(df))  # 1×N matrix of strings
  data = Matrix(df)
  writedlm(path, vcat(header, data), sep)
end

###################################################################################################

function writedf_dict(dir::String, dict::Dict{Int,DataFrame}; sep = ',')
  isdir(dir) || mkpath(dir)
  for (k, df) in dict
    path = joinpath(dir, "exp_$k.csv")
    writedf(path, df; sep = sep)
  end
end

###################################################################################################
