###################################################################################################
# HMGCR_TAM Metadata Extraction
###################################################################################################

begin
  # Load path definitions
  include(joinpath(PROGRAM_FILE === nothing ? "src" : "..", "config", "paths.jl"))
  using .Paths
  Paths.ensure_dirs()

  # Load configuration structs
  include(joinpath(Paths.CONFIG, "vars.jl"))
  include(joinpath(Paths.UTIL, "ioDataFrame.jl"))
  include(joinpath(Paths.UTIL, "ioLoadXLSX.jl"))
end;

###################################################################################################

# load packages
begin
    using DataFrames
    using XLSX
end;

###################################################################################################

"""
    extract_groups(path::String)

Open an XLSX file, read the "Metadata" sheet, and extract only the Group_* values.
Returns a DataFrame with columns: SourceFile, AnimalIndex, Group.
"""
function extract_groups(path::String)
    xf = XLSX.readxlsx(path)
    sheet = xf["Metadata"]

    # Flatten sheet into vector of strings
    rows = string.(vec(sheet[:]))

    # Find block delimiters
    start_idx = findfirst(contains("ANIMAL DATA REMARKS:"), rows)
    end_idx   = findfirst(contains("USER REMARKS:"), rows)

    if isnothing(start_idx) || isnothing(end_idx)
        @warn "No metadata block found in $path"
        return DataFrame(SourceFile=String[], AnimalIndex=Int[], Group=String[])
    end

    block = rows[start_idx+1:end_idx-1]

    groups = Union{String,Missing}[]
    for line in block
        tokens = split(line)
        # look for "Group_*" keys and take the following token as value
        for (i, tok) in enumerate(tokens)
            if startswith(tok, "Group_") && i < length(tokens)
                push!(groups, tokens[i+1])
            end
        end
    end

    # Pad to 16 animals
    while length(groups) < 16
        push!(groups, missing)
    end

    return DataFrame(SourceFile=fill(basename(path), 16),
                     AnimalIndex=1:16,
                     Group=groups)
end

###################################################################################################
# Main driver
###################################################################################################

files = filter(f -> endswith(f, ".xlsx"), readdir(Paths.HMGCR; join=true))

all_animals = vcat([extract_groups(f) for f in files]...)

# Summarize counts per group per file
group_summary = combine(groupby(all_animals, [:SourceFile, :Group]), nrow => :Count)

###################################################################################################
# Write outputs
###################################################################################################

# Write as tab‑separated files (or change sep="," for CSV)
writedf(joinpath(Paths.HMGCR, "all_animals.csv"), all_animals; sep = ',')
writedf(joinpath(Paths.HMGCR, "group_summary.csv"), group_summary; sep = ',')

###################################################################################################
