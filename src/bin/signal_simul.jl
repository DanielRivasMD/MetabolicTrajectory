using Random
using UnicodePlots

# Struct: only holds ranges and sampling setup
struct SignalParams
    duration::Float64                  # total seconds
    step::Float64                      # time step (1/fs)
    freq_range::Tuple{Float64,Float64} # (min Hz, max Hz)
    amp_range::Tuple{Float64,Float64}  # (min, max amplitude)
    phase_range::Tuple{Float64,Float64}# (min, max phase in radians)
end

# Generate one sine wave
function sine_wave(amplitude::Float64, frequency::Float64, phase::Float64, duration::Float64, step::Float64)
    t = 0:step:duration
    ω = 2π * frequency
    y = amplitude .* sin.(ω .* t .+ phase)
    return t, y
end

# Build composite signal: seed + num_components are external
function build_signal(sp::SignalParams, seed::Int, num_components::Int)
    Random.seed!(seed)
    t = 0:sp.step:sp.duration
    y_total = zeros(Float64, length(t))

    components = []
    for _ in 1:num_components
        f = rand() * (sp.freq_range[2] - sp.freq_range[1]) + sp.freq_range[1]
        a = rand() * (sp.amp_range[2] - sp.amp_range[1]) + sp.amp_range[1]
        p = rand() * (sp.phase_range[2] - sp.phase_range[1]) + sp.phase_range[1]
        push!(components, (frequency=f, amplitude=a, phase=p))
        _, y = sine_wave(a, f, p, sp.duration, sp.step)
        y_total .+= y
    end

    return t, y_total, components
end

# Example setup: 1000 samples at fs=1000 Hz
fs = 1000.0
n  = 1000
duration = n/fs
step = 1/fs

sp = SignalParams(duration, step, (10.0, 100.0), (1.0, 5.0), (0.0, 2π))

# Two reproducible signals with external seed + num_components
t1, y1, comps1 = build_signal(sp, 1234, 5)
t2, y2, comps2 = build_signal(sp, 5678, 5)

# Plot
lineplot(t1, y1; width=150, title="Signal 1 (seed=1234)") |> display
lineplot(t2, y2; width=150, title="Signal 2 (seed=5678)") |> display

# Inspect chosen components
println("Signal 1 components:")
for c in comps1
    println("freq=$(round(c.frequency, digits=2)) Hz, amp=$(round(c.amplitude, digits=2)), phase=$(round(c.phase, digits=2))")
end

println("\nSignal 2 components:")
for c in comps2
    println("freq=$(round(c.frequency, digits=2)) Hz, amp=$(round(c.amplitude, digits=2)), phase=$(round(c.phase, digits=2))")
end

