let

  # # interactive
  # include("src/live/saturatedSubsampling/config.jl")
  # include("src/live/saturatedSubsampling/util.jl")

  # Load path definitions
  include(joinpath((pwd() == @__DIR__) ? "src" : "../..", "config", "paths.jl"))
  using .Paths
  Paths.ensure_dirs()

  # Load configuration structs
  include(joinpath(Paths.CONFIG, "vars.jl"))
  include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
  include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))

  using Clustering
  using Dates
  using DataFrames
  using Distances
  using DynamicAxisWarping
  using LinearAlgebra
  using Plots
  using Statistics
  using Random
  using UMAP

  # Sigma experiment: one metadata XLSX, three CSV batches
  sigma_params = TrajectoryParams(
    metadata = Vars.SMETA_xlsx,
    batches = [Vars.SIG1R_HT_csv, Vars.SIG1R_WT_csv, Vars.KO_WT_csv],
  )

  # Load and split
  bundles = load_experiments(sigma_params)

  # Collect all metadata into one DataFrame
  meta = vcat([b.metadata for b in values(bundles)]...)

  # Keep only the columns of interest
  meta = select(meta, [:Animal_nr, :Sex, :Genotype])

  # Drop rows where Animal_nr == 0 (if those are placeholders)
  filter!(row -> row.Animal_nr != 0, meta)
  rename!(meta, :Animal_nr => :Animal)
  meta.Group = string.(meta.Sex, "_", meta.Genotype)
  animal_to_group = Dict(row.Animal => row.Group for row in eachrow(meta))

  subdfs = split_by_animal(bundles)

  sigma_test = [7, 8, 21, 22]
  sigma_dfs = Dict(k => subdfs[k] for k in sigma_test)

  using Plots
  using Colors

  function plot_RT_RER_one_by_one(subdfs)
    animals = sort(collect(keys(subdfs)))

    # compute global max length for x-axis
    maxlen = maximum(nrow(df) for df in values(subdfs))

    for animal in animals
      df = subdfs[animal]

      # metadata lookup
      row = meta[meta.Animal.==animal, :]
      sex = row.Sex[1]
      genotype = row.Genotype[1]

      # line color by sex
      linecolor = sex == "F" ? RGB(1, 0, 0) : RGB(0, 0, 1)

      # title color by genotype
      titlecolor =
        genotype == "S1RKO" ? RGB(0, 0, 0) :
        genotype == "Hetero" ? RGB(0.5, 0.5, 0.5) :
        genotype == "WT" ? RGB(0, 1, 0) : RGB(0.3, 0.3, 0.3)

      # vertical lines every 320 samples
      vlines = 320:320:maxlen

      plt = plot(
        df[!, "RT_RER"];
        color = linecolor,
        lw = 1.5,
        xlabel = "Time",
        ylabel = "RT_RER",
        ylim = (0, 2),
        xlim = (0, maxlen),
        title = "Animal $animal ($sex, $genotype)",
        titlefontcolor = titlecolor,
        legend = false,
        size = (800, 300),
      )

      # horizontal line at y = 0.7
      hline!(plt, [0.7], color = :black, alpha = 0.5, lw = 1)

      # vertical lines at multiples of 320
      vline!(plt, vlines, color = :black, alpha = 0.5, lw = 1)

      display(plt)
    end
  end

  plot_RT_RER_one_by_one(subdfs)

end
