# Chameleon

## Overview

`chameleon` contains a set of functions for automatically assigning colors to multi-dimensional data. The colors are
distinct (as much as possible) and also reflect the internal structure of the data (as much as possible). The function
for generating an arbitrary number of maximally-distinct colors is also useful by itself for generating color palettes
for arbitrary categorical data.

This is a port of the [published R package](https://CRAN.R-project.org/package=chameleon). Note the results aren't only
"mostly identical" for the same data, because the color conversion functions behave slightly differently between R and
Julia.

See the [v0.1.0 documentation](https://tanaylab.github.io/Chameleon.jl/v0.1.0) for details.

## Usage

In general, if you have a matrix whose rows represent some elements/observations and whose columns represent some
variables/measurements, then use `Chameleon.data_colors` to obtain a color for each row. You can pass this function the
optional `groups` parameter which assigns a group to each row, to obtain a color for each group instead of a color for
each row.

```julia
using Chameleon

data_matrix::AbstractMatrix{<:Real} = ...  # Arbitrary data

color_per_column::AbstractVector{<:AbstractString} =  # In #RRGGBB format
    data_colors(data_matrix)

group_per_column::AbstractVector{T} = ...  # Arbitrary type for group identifiers

color_per_group::AbstractDict{T, <:AbstractString} =  # Key is group identifier, value is in #RRGGBB format
    data_colors(data_matrix; groups = group_per_column)
```

If you just want a list of an arbitrary number of maximally-distinct colors, call `Chameleon.distinct_colors`. This will
return color names and their Lab coordinates, in an arbitrary order. You can use this to generate color palettes for
data without trying to make more-similar data have more-similar colors.

```julia
using Chameleon

n_colors::Integer = ...  # Non-negative

colors::AbstractVector{<:AbstractString}, # n_colors names in #RRGGBB format
labs::AbstractMatrix{<:AbstractFloat} = # 3 rows (l, a, b) and n_colors columns
    distinct_colors(n_colors)
```

## Installation

Just `Pkg.add("Chameleon")`, like installing any other Julia package.

## License (MIT)

Copyright © 2025 Weizmann Institute of Science

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
