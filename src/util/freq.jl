using UnicodePlots

# Define a struct to hold sine parameters
struct SineParams
  amplitude::Float64
  frequency::Float64   # Hz
  phase::Float64       # radians
  duration::Float64    # seconds
  step::Float64        # time step
end

# Function to generate the sine array
function sine_wave(params::SineParams)
  t = 0:params.step:params.duration
  ω = 2π * params.frequency
  y = params.amplitude .* sin.(ω .* t .+ params.phase)
  return t, y
end

# Example usage
params = SineParams(1.0, 1.0, 0.0, 10.0, 0.01)
t, y = sine_wave(params)

# Plot in terminal
lineplot(
  t,
  y;
  xlabel = "Time (s)",
  ylabel = "Amplitude",
  title = "Sine Function (UnicodePlots)",
)

