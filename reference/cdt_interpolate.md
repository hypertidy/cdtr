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

## Examples

``` r
P <- cbind(c(0, 1, 0), c(0, 0, 1)); T <- matrix(1:3, 1)
cdt_interpolate(P, T, cbind(v = c(10, 20, 30)), c(0.25, 0.5), c(0.25, 0.25))
#>         v
#> [1,] 17.5
#> [2,] 20.0
```
