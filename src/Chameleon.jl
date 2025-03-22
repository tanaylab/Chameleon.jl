"""
Functions for automatically assigning colors to multi-dimensional data. The colors are distinct (as much as possible)
and also reflect the internal structure of the data (as much as possible).
"""
module Chameleon

export data_colors
export distinct_colors

using Colors
using Distances
using Hungarian
using MultivariateStats
using StatsBase
using UMAP

"""
    data_colors(
        data::AbstractMatrix{<:Real};
        run_umap::Bool = true,
        minimal_saturation::Real = 33,
        minimal_lightness::Real = 20,
        maximal_lightness::Real = 80
    )::AbstractVector{<:AbstractString}

    data_colors(
        data::AbstractMatrix{<:Real};
        run_umap::Bool = true,
        groups::Union{AbstractVector{T}},
        minimal_saturation::Real = 33,
        minimal_lightness::Real = 20,
        maximal_lightness::Real = 80
    )::Dict{T, <:AbstractString}

Compute colors for multi-dimensional data.

Given a matrix of observation/element columns and variable/measurement rows, compute a color for each column (or group
of columns) such that the colors are distinct, and where more-similar colors roughly designate more-similar data columns
(or groups of columns).

This is intended to provide a "reasonable" set of colors to "arbitrary" data, for use as a convenient default when
investigating unknown data sets. It is not meant to replace hand-picked colors tailored for specific data (e.g., using
red colors for "bad" columns and green colors for "good" columns).

This ensures all colors are distinct by packing the (visible part) of the CIELAB color space with the needed number of
spheres using [`distinct_colors`](@ref). To assign the colors to the data, it uses UMAP to reduce the data to 3D. It
then uses principal component analysis to represent both the chosen colors (3D sphere centers) and the (3D UMAP) data as
point clouds with coordinates in the range 0-1, and finally uses a stable matching algorithm to map these point clouds
to each other, thereby assigning a color to each data column. If the data is grouped, then the center of gravity of each
group is used to generate a color for each group.

# Arguments

  - `data`: A matrix whose columns represent elements/observations and columns represent
    variables/measurements.

  - `run_umap`: A Boolean specifying whether to run UMAP on the data to convert it to 3D. If `false`, the data matrix must
    have exactly 3 columns and will be used as-is.
  - `groups`: An optional array with an entry per column containing the identifier of the group the column belongs to.
  - `minimal_saturation`: Exclude colors whose saturation (hypot(a, b) in CIELAB color space) is less than this value.
  - `minimal_lightness`: Exclude colors whose lightness (l in CIELAB color space) is less than this value.
  - `maximal_lightness`: Exclude colors whose lightness (l in CIELAB color space) is more than this value.

# Returns

If `groups` was specified, an dictionary mapping the group names to colors. Otherwise, a vector with one entry per
column containing its color. Colors are named in the usual `#RRGGBB` hexadecimal format.

!!! note

    Since this calls UMAP, which (sometimes) uses randomness under the hood, and doesn't have a way to pass it a random
    number generator, then if you want reproducible results, you need to call `Random.seed!(...)` yourself immediately
    before calling `data_colors`. If/when UMAP allows passing the random number generator as a parameter (as it should),
    then we should add the same parameter here and pass it to UMAP.
"""
function data_colors(
    data::AbstractMatrix{<:Real};
    run_umap::Bool = true,
    groups::Union{AbstractVector, Nothing} = nothing,
    minimal_saturation::Real = 33,
    minimal_lightness::Real = 20,
    maximal_lightness::Real = 80,
)::Union{Dict{<:Any, <:AbstractString}, AbstractVector{<:AbstractString}}
    n_rows, n_columns = size(data)
    if groups !== nothing
        @assert length(groups) == n_columns
    end
    if !run_umap
        @assert n_rows == 3
    end

    if n_columns == 0
        if groups === nothing
            return AbstractString[]
        else
            return Dict{eltype(groups), AbstractString}()
        end
    end

    if groups !== nothing
        data, unique_groups = group_centers(data, groups)
        n_groups = length(unique_groups)
    else
        unique_groups = nothing
        n_groups = n_columns
    end

    colors = distinct_colors(n_groups; minimal_saturation, minimal_lightness, maximal_lightness)

    if n_groups <= 2
        result = colors[1]

    else
        @assert run_umap || size(data) == (3, n_groups)
        if run_umap && n_rows > 3
            n_neighbors = min(15, n_groups - 1)
            data = umap(data, 3; n_neighbors)  # NOJET
        end

        data_pca = fit(PCA, data; maxoutdim = 3)
        data_points = normalized_data(MultivariateStats.predict(data_pca, data))

        color_pca = fit(PCA, colors[2]; maxoutdim = 3)  # NOJET
        color_points = normalized_data(MultivariateStats.predict(color_pca, colors[2]))

        distances = pairwise(Euclidean(), data_points, color_points; dims = 2)
        @assert size(distances) == (n_groups, n_groups)

        assignment = hungarian(distances)[1]
        @assert length(assignment) == n_groups

        result = colors[1][assignment]
    end

    @assert length(result) == n_groups

    if unique_groups === nothing
        return result
    else
        return Dict(zip(unique_groups, result))
    end
end

function normalized_data(data::AbstractMatrix{<:AbstractFloat})::AbstractMatrix{<:AbstractFloat}
    n_rows, n_columns = size(data)
    if n_rows == 1
        data = vcat(data, collect(1:n_columns)', collect(n_columns:-1:1)')
    elseif size(data, 1) == 2
        data = vcat(data, collect(1:n_columns)')
    end
    @assert size(data) == (3, n_columns)

    low = vec(minimum(data; dims = 2))
    high = vec(minimum(data; dims = 2))
    range = high .- low
    @assert length(range) == 3

    range[range .== 0] .= 1
    result = (data .- low) ./ range
    @assert size(result) == (3, n_columns)

    return result
end

function group_centers(
    data::AbstractMatrix{<:Real},
    groups::AbstractVector{T},
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractVector{T}} where {T}
    unique_groups = unique(groups)

    n_groups = length(unique_groups)
    n_rows = size(data, 1)

    grouped_data = Matrix{Float64}(undef, n_rows, n_groups)

    for (group_index, group_value) in enumerate(unique(groups))
        group_mask = groups .== group_value
        @assert any(group_mask)
        grouped_data[:, group_index] .= mean(filter(!isnan, data[:, group_mask]))
    end

    return grouped_data, unique_groups
end

"""
    distinct_colors(
        n_colors::Integer;
        minimal_saturation::Real=33,
        minimal_lightness::Real=20,
        maximal_lightness::Real=80
    )::Tuple{AbstractVector{<:AbstractString}, AbstractMatrix{<:Real}}

Pick a number of distinct colors.

This ensures all colors are distinct by packing the (visible part) of the CIELAB color space with the needed number of
spheres, and using their centers to generate the colors.

# Arguments

  - `n_colors`: The requested (positive) number of colors.

  - `minimal_saturation`: Exclude colors whose saturation (hypot(a, b) in CIELAB color space) is less than this value.
  - `minimal_lightness`: Exclude colors whose lightness (l in CIELAB color space) is less than this value.
  - `maximal_lightness=80`: Exclude colors whose lightness (l in CIELAB color space) is more than this value.

# Returns

A tuple with two elements, a vector of color names and a matrix of the `lab` coordinates of the colors (rows are `l`,
`a` and `b`, columns are colors). Colors are named in the usual `#RRGGBB` hexadecimal format.
"""
function distinct_colors(
    n_colors::Integer;
    minimal_saturation::Real = 33,
    minimal_lightness::Real = 20,
    maximal_lightness::Real = 80,
)::Tuple{AbstractVector{<:AbstractString}, AbstractMatrix{<:Real}}
    @assert n_colors > 0
    @assert 0 < minimal_saturation < 150
    @assert 0 <= minimal_lightness <= maximal_lightness <= 100

    too_large_step = 100.0
    large_step = 100.0
    large_step_colors = pick_step_colors(large_step, minimal_saturation, minimal_lightness, maximal_lightness)

    while large_step_colors === nothing
        too_large_step = large_step
        large_step = large_step / 2.0
        large_step_colors = pick_step_colors(large_step, minimal_saturation, minimal_lightness, maximal_lightness)
    end

    if length(large_step_colors[1]) == n_colors
        return large_step_colors  # UNTESTED
    end

    if length(large_step_colors[1]) > n_colors
        while true
            if length(large_step_colors[1]) == 4
                large_step_colors_names = large_step_colors[1][1:n_colors]
                large_step_colors_labs = large_step_colors[2][:, 1:n_colors]
                return (large_step_colors_names, large_step_colors_labs)
            end

            mid_step = (too_large_step + large_step) / 2
            mid_step_colors = pick_step_colors(mid_step, minimal_saturation, minimal_saturation, maximal_lightness)

            if mid_step_colors === nothing
                too_large_step = mid_step
                continue
            end

            if length(mid_step_colors[1]) == n_colors
                return mid_step_colors
            end

            large_step = mid_step
            large_step_colors = mid_step_colors
        end
    end

    @assert length(large_step_colors[1]) < n_colors

    small_step = large_step
    small_step_colors = large_step_colors
    while length(small_step_colors[1]) < n_colors
        small_step = small_step / 2
        small_step_colors = pick_step_colors(small_step, minimal_saturation, minimal_lightness, maximal_lightness)
        @assert !isnothing(small_step_colors)
    end

    @assert length(large_step_colors[1]) < n_colors
    @assert length(small_step_colors[1]) >= n_colors

    while length(small_step_colors[1]) > n_colors
        if large_step - small_step < 1e-6
            small_step_colors_names = small_step_colors[1][1:n_colors]
            small_step_colors_labs = small_step_colors[2][:, 1:n_colors]
            small_step_colors = (small_step_colors_names, small_step_colors_labs)
            break
        end

        mid_step = (small_step + large_step) / 2
        mid_step_colors = pick_step_colors(mid_step, minimal_saturation, minimal_lightness, maximal_lightness)
        @assert mid_step_colors !== nothing
        if length(mid_step_colors[1]) >= n_colors
            small_step = mid_step
            small_step_colors = mid_step_colors
        else
            large_step = mid_step
            large_step_colors = mid_step_colors
        end
    end

    @assert length(small_step_colors[1]) == n_colors
    return small_step_colors
end

function pick_step_colors(
    step::AbstractFloat,
    minimal_saturation::Real,
    minimal_lightness::Real,
    maximal_lightness::Real,
)::Union{Tuple{AbstractVector{<:AbstractString}, AbstractMatrix{<:AbstractFloat}}, Nothing}
    lab = lab_tetragrid(step)

    saturation = sqrt.(lab[2, :] .^ 2 .+ lab[3, :] .^ 2)
    mask = (saturation .>= minimal_saturation) .& (lab[1, :] .>= minimal_lightness) .& (lab[1, :] .<= maximal_lightness)  # NOJET
    if sum(mask) < 4
        return nothing  # UNTESTED
    end

    lab = lab[:, mask]

    srgb = map(column -> safe_convert_lab_to_srgb(column), eachcol(lab))  # NOJET
    mask = srgb .!== nothing
    if sum(mask) < 4
        return nothing
    end

    lab = lab[:, mask]
    srgb = srgb[mask]

    colors = map(color -> "#" * hex(color), srgb)

    return (colors, lab)
end

function safe_convert_lab_to_srgb(lab::AbstractVector{<:AbstractFloat})::Union{RGB, Nothing}
    srgb = convert(RGB, Lab(lab...))
    back_lab = convert(Lab, srgb)
    if maximum(abs.(lab .- [back_lab.l, back_lab.a, back_lab.b])) < 1e-6
        return srgb
    else
        return nothing
    end
end

function lab_tetragrid(step::Real)::AbstractMatrix{<:AbstractFloat}
    l_steps = round(Int, 100 / (step * 2 * sqrt(6) / 3))
    a_steps = round(Int, 250 / (step * 2))
    b_steps = round(Int, 250 / (step * sqrt(3)))
    grid = Matrix{Float64}(undef, 3, (l_steps + 1) * (a_steps + 1) * (b_steps + 1))

    i = 1
    for li in 0:l_steps
        for la in 0:a_steps
            for lb in 0:b_steps
                grid[1, i] = step / 2 + (li * 2 * sqrt(6) / 3) * step
                grid[2, i] = step / 2 - 100 + (2 * la + (lb + li) % 2) * step
                grid[3, i] = step / 2 - 150 + (sqrt(3) * (lb + (li % 2) / 3)) * step
                i += 1
            end
        end
    end

    return grid
end

end # module
