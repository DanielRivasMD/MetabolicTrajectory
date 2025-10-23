let

  using DynamicAxisWarping
  using Plots
  using Random

  # 1) Seed RNG for reproducibility
  Random.seed!(1234)

  # 2) Generate two random sequences
  len = 200
  seq1 = randn(len)                  # random Gaussian noise
  seq2 = randn(len) .+ 0.5           # shifted random noise

  # Plot the raw sequences
  plt_seq = plot(seq1; label = "seq1", color = :blue, lw = 2)
  plot!(plt_seq, seq2; label = "seq2", color = :red, lw = 2, title = "Random Sequences")
  display(plt_seq)

  # 3) Compute DTW cost matrix
  D_rand = dtw_cost_matrix(seq1, seq2, SqEuclidean())

  # 4a) Heatmap of cost matrix
  plt_heat_rand = heatmap(
    D_rand;
    xlabel = "seq1 index",
    ylabel = "seq2 index",
    title = "DTW Cost Matrix (Random Sequences)",
    colorbar_title = "Cost",
  )
  display(plt_heat_rand)

  # 4b) Alignment overlay on original signals
  cost_rand, i1_rand, i2_rand = dtw(seq1, seq2, SqEuclidean())

  plt_align_rand = plot(seq1; label = "seq1", color = :blue, lw = 2)
  plot!(plt_align_rand, seq2; label = "seq2", color = :red, lw = 2)

  # Draw alignment lines (subsampled for clarity)
  for (idx1, idx2) in zip(i1_rand[1:10:end], i2_rand[1:10:end])
    plot!(
      plt_align_rand,
      [idx1, idx2],
      [seq1[idx1], seq2[idx2]];
      color = :gray,
      alpha = 0.4,
      label = false,
    )
  end

  display(plt_align_rand)

  println("Euclidean distance (seq1 vs seq2): ", evaluate(Euclidean(), seq1, seq2))
  println("DTW cost (seq1 vs seq2): ", cost_rand)

end
