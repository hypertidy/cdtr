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
cdt_locate(P, T, xq, yq)
```

## Arguments

- P:

  vertex coordinates (n x 2)

- T:

  triangle vertex indices (m x 3, 1-based)

- xq, yq:

  query coordinates

## Value

an integer vector, `length(xq)`, of 1-based row indices into `T` (`NA`
outside the mesh)

## Details

A point on an edge or vertex shared by several triangles is given to the
lowest triangle index that contains it. Points outside every triangle,
and `NA` points, give `NA`. Triangles may be in either orientation, and
degenerate (zero area) triangles contain no point.

## Examples

``` r
P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
T <- rbind(c(1, 2, 3), c(1, 3, 4))
## (0.5, 0.5) is on the shared diagonal: the lower triangle id wins
cdt_locate(P, T, xq = c(0.75, 0.25, 0.5, 2), yq = c(0.25, 0.75, 0.5, 2))
#> [1]  1  2  1 NA
```
