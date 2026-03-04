####################################################################################################
# src/core.jl
####################################################################################################

module Process

####################################################################################################

using Dates, DataFrames, Distances, DynamicAxisWarping, Statistics

####################################################################################################

struct SubSampleID
  subject::Int
  ixs::Tuple{Int,Int}
  time::Tuple{DateTime,DateTime}
end

struct SubSampleContainer
  subsamples::Vector{Vector{Float64}}
  ids::Vector{SubSampleID}
end

struct TrajectoryParams
  metadata::String
  batches::Vector{String}
  meta_path::String
  dtw_path::String
  nsamples::Int
  len::Float64
  var::Float64
  limits::Tuple{Int,Int}
end

####################################################################################################

diff_cumulative(v::Vector{<:Real}) = [v[1]; v[2:end] .- v[1:end-1]]

# function collect_subsamples(
#   signal::Vector{Float64},
#   times::Vector{DateTime},
#   params::TrajectoryParams,
#   diff::Bool,
# )
#   # ... your existing logic ...
#   return SubSampleContainer(subsamples, ids)
# end

####################################################################################################

function compute_dtw_matrix(signals::Vector{Vector{Float64}})
  N = length(signals)
  cost = zeros(N, N)
  for i = 1:N, j = i:N
    d, _, _ = dtw(signals[i], signals[j], SqEuclidean())
    nd = d / mean([length(signals[i]), length(signals[j])])
    cost[i, j] = cost[j, i] = nd
  end
  return cost
end

####################################################################################################

function readdf(path; sep = '\t')
  data, header = readdlm(path, sep, header = true)
  return DataFrame(data, vec(header))
end

function writedf(path, df::DataFrame; sep = ',')
  header = permutedims(names(df))
  data = Matrix(df)
  writedlm(path, vcat(header, data), sep)
end

####################################################################################################

end

####################################################################################################

