####################################################################################################

using ArgParse

####################################################################################################

const HELP =
  "\e[1;32mDaniel Rivas\e[0m " * "\e[3;90m<danielrivasmd@gmail.com>\e[0m\n\n\n\n\n"

####################################################################################################

"Convert struct fields to a Dict with Symbol keys"
function struct_to_dict(x)
  Dict(name => getfield(x, name) for name in propertynames(x))
end

"Convert Dict with String keys (from TOML) into Dict with Symbol keys"
function symbolise_keys(d::Dict)
  Dict(Symbol(k) => v for (k, v) in d)
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
