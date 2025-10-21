let

# 1) Seed RNG
Random.seed!(2025)

# 2) Generate noisy sinusoidal signals
x = range(0, 4π; length=400)
signal1 = sin.(x) .+ 0.1 .* randn(length(x))   # sine + small noise
signal2 = cos.(x) .+ 0.1 .* randn(length(x))   # cosine + small noise

# Plot raw noisy signals
plt_signals = plot(x, signal1; label="noisy sin(x)", color=:blue, lw=2)
plot!(plt_signals, x, signal2; label="noisy cos(x)", color=:red, lw=2,
      title="Noisy Sinusoidal Signals")
display(plt_signals)

# 3) Compute DTW cost matrix
D_noisy = dtw_cost_matrix(signal1, signal2, SqEuclidean())

# 4a) Heatmap of cost matrix
plt_heat_noisy = heatmap(D_noisy;
    xlabel="sin(x) index",
    ylabel="cos(x) index",
    title="DTW Cost Matrix (Noisy Sinusoids)",
    colorbar_title="Cost"
)
display(plt_heat_noisy)

# 4b) Alignment overlay on original signals
cost_noisy, i1_noisy, i2_noisy = dtw(signal1, signal2, SqEuclidean())

plt_align_noisy = plot(x, signal1; label="noisy sin(x)", color=:blue, lw=2)
plot!(plt_align_noisy, x, signal2; label="noisy cos(x)", color=:red, lw=2)

# Draw alignment lines (subsampled for clarity)
for (idx1, idx2) in zip(i1_noisy[1:15:end], i2_noisy[1:15:end])
    plot!(plt_align_noisy, [x[idx1], x[idx2]], [signal1[idx1], signal2[idx2]];
          color=:gray, alpha=0.4, label=false)
end

display(plt_align_noisy)

println("Euclidean distance (noisy sin vs noisy cos): ",
        evaluate(Euclidean(), signal1, signal2))
println("DTW cost (noisy sin vs noisy cos): ", cost_noisy)

end
