###################################################################################################

using DataFrames
using Dates
using DelimitedFiles
using XLSX

##################################################################################################

function castdf!(df::DataFrame)
  # 1) first column → DateTime
  df[!, 1] = DateTime.(df[!, 1])

  # 2) all other columns → Float64
  for c in names(df)[2:end]
    df[!, c] = Float64.(df[!, c])
  end

  return df
end

##################################################################################################

function readdf(path; sep = '\t')
  data, header = readdlm(path, sep, header = true)
  return DataFrame(data, vec(header))
end

###################################################################################################

function writedf(path, df::DataFrame; sep = ',')
  header = permutedims(names(df))  # 1×N matrix of strings
  data = Matrix(df)
  writedlm(path, vcat(header, data), sep)
end

###################################################################################################

function readdf_dict(dir::String; sep = ',')
  dict = Dict{Int,DataFrame}()

  for file in readdir(dir)
    m = match(r"^exp_(\d+)\.csv$", file)
    m === nothing && continue

    key = parse(Int, m.captures[1])
    df = readdf(joinpath(dir, file); sep = sep)
    castdf!(df)
    dict[key] = df
  end

  return dict
end

###################################################################################################

function writedf_dict(dir::String, dict::Dict{Int,DataFrame}; sep = ',')
  for (k, df) in dict
    path = joinpath(dir, "exp_$k.csv")
    writedf(path, df; sep = sep)
  end
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
