####################################################################################################

using DynamicAxisWarping

####################################################################################################

function dtw_distance_matrix_by_variable(df::DataFrame, var::Symbol)
    gdf = groupby(df, :Animal)

    # keep animal IDs as strings for labeling
    series = Dict(string(first(subdf.Animal)) => collect(skipmissing(subdf[!, var]))
                  for subdf in gdf)

    animals = sort(collect(keys(series)))
    n = length(animals)
    M = fill(NaN, n, n)

    for i in 1:n, j in i:n
        si, sj = series[animals[i]], series[animals[j]]
        isempty(si) || isempty(sj) && continue
        d = dtw(si, sj)
        M[i,j] = d
        M[j,i] = d
    end

    return animals, M
end

####################################################################################################

using DataFrames

function dtw_distance_matrix_by_variable(df::DataFrame, var::Symbol;
                                         dist=(a,b)->abs(a-b),
                                         window::Union{Nothing,Int}=nothing,
                                         normalize::Bool=true)
    gdf = groupby(df, :Animal)
    animals = String.(unique(df.Animal))
    n = length(animals)

    # map animal -> series
    series = Dict{String,Vector{Float64}}()
    for subdf in gdf
        a = string(first(subdf.Animal))
        s = series_from(sort_by_time(subdf), var)
        isempty(s) && continue
        series[a] = s
    end

    animals = sort(collect(keys(series)))
    n = length(animals)
    M = fill(NaN, n, n)

    for i in 1:n, j in i:n
        si = series[animals[i]]
        sj = series[animals[j]]
        if isempty(si) || isempty(sj)
            M[i,j] = NaN; M[j,i] = NaN
            continue
        end
        d, _ = dtw(si, sj; dist=dist, window=window, normalize=normalize)
        M[i,j] = d
        M[j,i] = d
    end

    return animals, M
end

####################################################################################################

# Minimal, robust DTW (L1/L2 selectable) with path reconstruction
function dtw(s::AbstractVector{<:Real}, t::AbstractVector{<:Real};
             dist::Function = (a,b)->abs(a-b),  # L1 default
             window::Union{Nothing,Int}=nothing,
             normalize::Bool=true)
    ns, nt = length(s), length(t)
    D = fill(Inf, ns+1, nt+1)
    D[1,1] = 0.0

    # optional Sakoe-Chiba band
    if window !== nothing
        w = max(window, abs(ns - nt))
    else
        w = max(ns, nt) # effectively no constraint
    end

    for i in 1:ns
        jmin = max(1, i - w)
        jmax = min(nt, i + w)
        @inbounds for j in jmin:jmax
            cost = dist(s[i], t[j])
            D[i+1, j+1] = cost + min(D[i, j+1], D[i+1, j], D[i, j]) # diagonal, left, up
        end
    end

    # backtrack path
    i, j = ns+1, nt+1
    path = Vector{Tuple{Int,Int}}()
    while i>1 && j>1
        push!(path, (i-1, j-1))
        # choose predecessor
        a, b, c = D[i-1, j], D[i, j-1], D[i-1, j-1]
        if c <= a && c <= b
            i -= 1; j -= 1
        elseif a <= b
            i -= 1
        else
            j -= 1
        end
    end
    while i>1; push!(path, (i-1, j-1)); i -= 1; end
    while j>1; push!(path, (i-1, j-1)); j -= 1; end
    reverse!(path)

    distval = D[ns+1, nt+1]
    if normalize
        distval /= (ns + nt)  # length‑normalized distance
    end
    return distval, path
end

####################################################################################################

# Extract numeric vector from a grouped subdf for a given variable
function series_from(subdf::DataFrame, var::Symbol)
    v = subdf[!, var]
    v = collect(skipmissing(v))
    return Float64.(v)
end

# Optional: restrict by DateTime window or sort by time if needed
function sort_by_time(subdf::DataFrame)
    sort(subdf, :DateTime)
end

####################################################################################################

using UnicodePlots

function plot_dtw_matrix_unicode(animals::Vector{String}, M::Matrix{Float64}, var::Symbol)
    # Convert to a simple image plot in the terminal
    display(heatmap(1:length(animals), 1:length(animals), M;
                    xlabel="Animal index",
                    ylabel="Animal index",
                    title="DTW distances for $(string(var))"))
    # Also print a tidy table
    println("Animals order: ", animals)
end

function plot_alignment_unicode(df::DataFrame, var::Symbol, a1::AbstractString, a2::AbstractString;
                                window::Union{Nothing,Int}=nothing)
    s1 = series_from(sort_by_time(df[df.Animal .== a1, :]), var)
    s2 = series_from(sort_by_time(df[df.Animal .== a2, :]), var)
    d, path = dtw(s1, s2; window=window)

    # Build aligned sequences using path
    y1 = Float64[]; y2 = Float64[]
    for (i,j) in path
        push!(y1, s1[i]); push!(y2, s2[j])
    end

    # Plot both aligned sequences against path index
    idx = 1:length(y1)
    plt = lineplot(idx, y1; name="$(a1)", xlabel="Warp index", ylabel=string(var),
                   title="DTW align $(a1) vs $(a2) | d=$(round(d, digits=4))")
    lineplot!(plt, idx, y2; name="$(a2)")
    display(plt)
end

####################################################################################################

using Statistics

function zscore(v::Vector{Float64})
    μ = mean(v); σ = std(v)
    σ == 0 ? fill(0.0, length(v)) : (v .- μ) ./ σ
end

function series_from_scaled(subdf::DataFrame, var::Symbol)
    v = Float64.(collect(skipmissing(subdf[!, var])))
    return zscore(v)
end

####################################################################################################
