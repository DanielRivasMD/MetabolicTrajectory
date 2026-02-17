####################################################################################################

using ArgParse

####################################################################################################

const HELP =
  "\e[1;32mDaniel Rivas\e[0m " * "\e[3;90m<danielrivasmd@gmail.com>\e[0m\n\n\n\n\n"

####################################################################################################

function struct_to_dict(x)
  Dict(name => getfield(x, name) for name in propertynames(x))
end

function symbolise_keys(d::Dict)
  Dict{Symbol,Any}(Symbol(k) => v for (k, v) in d)
end

####################################################################################################

function sigma_args()
  desc = HELP * "Extract a random subsample & calculate Dynamic Time Wraping matrix\n"
  s = ArgParseSettings(description = desc)

  @add_arg_table s begin
    "--config"
    help = "Path to config TOML"

    "--verbose", "-v"
    help = "Enable verbose logging"
    action = :store_true
  end

  return parse_args(s)
end

####################################################################################################
