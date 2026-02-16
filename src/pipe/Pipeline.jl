module Pipeline

export Stage, run_pipeline

using Logging

struct Stage
  name::String
  func::Function
end

"""
    run_pipeline(stages; params=Dict())

Runs each stage in order, logging start/end times.
"""
function run_pipeline(stages; params = Dict())
  @info "Starting pipeline with $(length(stages)) stages"

  for s in stages
    @info "→ Running stage: $(s.name)"
  try
      s.func(params)
      @info "✓ Finished stage: $(s.name)"
    catch e
      @error "Stage $(s.name) failed" exception = e
      rethrow(e)
    end
  end

  @info "Pipeline completed"
end

end
