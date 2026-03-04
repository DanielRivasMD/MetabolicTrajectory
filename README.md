MetabolicTrajetories/
├── Project.toml          # depends on Avicenna, plus any other packages
├── src/
│   ├── core.jl           # Pure analysis logic (types, functions) – the “what”
│   ├── workflows.jl      # Workflow definitions using Avicenna.Workflow – the “how”
│   ├── cli.jl            # CLI entry points using Avicenna.CLI
│   ├── repl.jl           # Interactive helpers using Avicenna.REPL
│   ├── document.jl       # Report generation using Avicenna.Document
│   └── utils/
│       ├── io.jl         # I/O helpers (readdf, writedf, etc.)
│       ├── params.jl     # Parameter loading (TrajectoryParams, etc.)
│       └── wrangle.jl    # Data wrangling (split_by_animal, etc.)
├── bin/
│   ├── run_sigma.jl      # CLI entry (calls your CLI module)
│   └── run_plot.jl       # Separate plotting CLI
├── scripts/              # (optional) legacy wrappers can call the new CLI
├── test/
└── README.md
