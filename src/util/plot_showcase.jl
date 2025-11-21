using Plots

# Dummy data
data = rand(100, 100)

# List of palettes to showcase
palettes = [
  :viridis,
  :plasma,
  :inferno,
  :magma,
  :cividis,
  :thermal,
  :balance,
  :Spectral,
  :Dark2_8,
  :tab10,
]

# Create one heatmap per palette
plots = [
  heatmap(data; color = cgrad(p, 200), title = string(p), legend = false, ticks = false)
  for p in palettes
]

# Arrange in a 2×5 grid
plot(plots...; layout = (2, 5), size = (2000, 800))

