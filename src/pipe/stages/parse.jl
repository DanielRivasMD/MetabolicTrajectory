module ParseStage

export run_parse

using Logging

function run_parse(params)
  @info "[parse] reading raw data"
  sleep(0.5)  # mock work
  @info "[parse] done"
end

end
