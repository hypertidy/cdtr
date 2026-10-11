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
## two triangles on the unit square, z = x + y at the corners
P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
T <- rbind(c(1, 2, 3), c(1, 3, 4))
A <- cbind(z = P[, 1] + P[, 2])
## a linear field is reproduced exactly inside; (2, 2) is outside
cdt_interpolate(P, T, A, xq = c(0.5, 0.25, 2), yq = c(0.5, 0.1, 2))
#>         z
#> [1,] 1.00
#> [2,] 0.35
#> [3,] 2.00
```
