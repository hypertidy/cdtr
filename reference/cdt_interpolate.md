# Interpolate per-vertex attributes onto points

Linear (barycentric) interpolation of a matrix of per-vertex values from
a triangle mesh onto query points. Points in no triangle take the
nearest vertex's value. Used to carry attributes such as z onto Steiner
vertices after refinement, matching the `PA` behaviour of RTriangle.

## Usage

``` r
cdt_interpolate(P, T, A, xq, yq)
```

## Arguments

- P:

  vertex coordinates (n x 2)

- T:

  triangle vertex indices (m x 3, 1-based)

- A:

  per-vertex attribute matrix (n x k) or vector

- xq, yq:

  query coordinates

## Value

a length(xq) x k matrix (column names kept from `A`)
