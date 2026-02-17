####################################################################################################

using ArgParse

####################################################################################################

const HELP =
  "\e[1;32mDaniel Rivas\e[0m " * "\e[3;90m<danielrivasmd@gmail.com>\e[0m\n\n\n\n\n"

####################################################################################################

macro vinfo(args, msg)
  return :($args["verbose"] ? @info($msg) : nothing)
end

function struct_to_dict(x)
  Dict(name => getfield(x, name) for name in propertynames(x))
end

function symbolise_keys(d::Dict)
  Dict{Symbol,Any}(Symbol(k) => v for (k, v) in d)
end

####################################################################################################

function sigma_args()
  desc = HELP * "Extract a random subsample & calculate Dynamic Time Wraping\n"
  s = ArgParseSettings(description = desc)

  @add_arg_table s begin
    "--trajectory"
    help = "Path to Trajectory TOML"

    "--verbose", "-v"
    help = "Enable verbose logging"
    action = :store_true
  end

  return parse_args(s)
end

####################################################################################################

function pipe_args()
  desc = HELP * "Pipeline"
  s = ArgParseSettings(description = desc)

  @add_arg_table s begin
    "--config"
    help = "Optional TOML config file for PipeParams"
    arg_type = String

    "--input"
    help = "Input directory"
    arg_type = String

    "--output"
    help = "Output directory"
    arg_type = String

    "--verbose", "-v"
    help = "Enable verbose logging"
    action = :store_true
  end

  return parse_args(s)
end

####################################################################################################
