####################################################################################################

using ArgParse

####################################################################################################

const HELP =
  "\e[1;32mDaniel Rivas\e[0m " * "\e[3;90m<danielrivasmd@gmail.com>\e[0m\n\n\n\n\n"

####################################################################################################

function struct_to_dict(x)
  d = Dict{Symbol,Any}()
  for name in propertynames(x)
    d[name] = getfield(x, name)
  end
  return d
end

function symbolise_keys(d::Dict)
  out = Dict{Symbol,Any}()
  for (k, v) in d
    out[Symbol(k)] = v
  end
  return out
end

####################################################################################################

function subsample_args()
  desc = HELP * "Extract a random subsample & calculate Dynamic Time Wraping\n"
  s = ArgParseSettings(description = desc)

  @add_arg_table s begin
    "--trajectory"
    help = "Path to Trajectory TOML"
  end

  return parse_args(s)
end

####################################################################################################
