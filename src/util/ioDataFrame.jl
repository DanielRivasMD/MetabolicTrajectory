###################################################################################################

# load packages
begin
  using DataFrames
  using DelimitedFiles
  using XLSX
end;

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

function split_by_animal(df::DataFrame; timecol::Symbol = :Date_Time)
  colnames = names(df)

  # Map each column to its animal ID
  animal_ids = Dict{Symbol,Int}()
  for col in colnames
    m = match(r"_(\d+)$", String(col))
    if m !== nothing
      animal_ids[Symbol(col)] = parse(Int, m.captures[1])
    end
  end

  # Group columns by animal ID
  grouped = Dict{Int,Vector{Symbol}}()
  for (col, id) in animal_ids
    push!(get!(grouped, id, Symbol[]), col)
  end

  subdfs = Dict{Int,DataFrame}()
  for (id, cols) in grouped
    sdf = df[:, vcat([timecol], sort(cols))]

    # Rename columns to strip suffix
    renames = Dict(c => Symbol(replace(String(c), r"_\d+$" => "")) for c in cols)
    rename!(sdf, renames)

    # Reorder: time column first, then sorted base names for this animal
    bases = sort(Symbol.(replace.(String.(cols), r"_\d+$" => "")))
    subdfs[id] = sdf[:, vcat([timecol], bases)]
  end

  return subdfs
end

###################################################################################################
