# Locate points in a triangle mesh

Find the triangle that contains each query point. This is the point
locator inside
[`cdt_interpolate()`](https://hypertidy.github.io/cdtr/reference/cdt_interpolate.md)
(a uniform grid of triangle buckets), with exact orientation tests from
the bundled robust predicates, so a point is never lost to rounding
between two triangles. With a triangle to feature table it answers
point-in-polygon for many points at once.

## Usage

``` r
cdt_locate(P, T, xq, yq, barycentric = FALSE)
```

## Arguments

- P:

  vertex coordinates (n x 2)

- T:

  triangle vertex indices (m x 3, 1-based)

- xq, yq:

  query coordinates

- barycentric:

  if `TRUE`, return a data frame with the barycentric coordinates,
  rather than only the triangle index

## Value

by default an integer vector, `length(xq)`, of 1-based row indices into
`T` (`NA` outside the mesh). With `barycentric = TRUE`, a data frame
with one row per query point and columns `triangle` (the same index),
and `b1`, `b2`, `b3`, which are `NA` where `triangle` is.

## Details

A point on an edge or vertex shared by several triangles is given to the
lowest triangle index that contains it. Points outside every triangle,
and `NA` points, give `NA`. Triangles may be in either orientation, and
degenerate (zero area) triangles contain no point.

With `barycentric = TRUE` the barycentric coordinates of each point in
its triangle are returned as well. They are the weights of the
triangle's three vertices (in the order given in `T`) that reproduce the
point, so `b1 * P[T[i, 1], ] + b2 * P[T[i, 2], ] + b3 * P[T[i, 3], ]` is
the point, and any per-vertex value is interpolated by the same weights,
as
[`cdt_interpolate()`](https://hypertidy.github.io/cdtr/reference/cdt_interpolate.md)
does. The weights are not negative and sum to 1; they are computed in
double precision, with rounding error clamped, so a point on an edge has
a weight of exactly 0 for the opposite vertex.

## Examples

``` r
P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
T <- rbind(c(1, 2, 3), c(1, 3, 4))
## (0.5, 0.5) is on the shared diagonal: the lower triangle id wins
cdt_locate(P, T, xq = c(0.75, 0.25, 0.5, 2), yq = c(0.25, 0.75, 0.5, 2))
#> [1]  1  2  1 NA

## the weights of each triangle's vertices
loc <- cdt_locate(P, T, xq = c(0.75, 0.25, 0.5), yq = c(0.25, 0.75, 0.5),
                  barycentric = TRUE)
loc
#>   triangle   b1   b2   b3
#> 1        1 0.25 0.50 0.25
#> 2        2 0.25 0.25 0.50
#> 3        1 0.50 0.00 0.50
## a value known at the vertices, interpolated with those weights
z <- c(0, 10, 20, 10)
rowSums(as.matrix(loc[c("b1", "b2", "b3")]) * matrix(z[T[loc$triangle, ]], ncol = 3))
#> [1] 10 10 10
```
