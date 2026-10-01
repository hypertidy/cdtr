# Unique edges of a triangulation

The index pairs of every triangle edge once, for plotting with
[`graphics::segments()`](https://rdrr.io/r/graphics/segments.html) or
for building an edge table.

## Usage

``` r
cdt_edges(x)
```

## Arguments

- x:

  result of
  [`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)

## Value

integer matrix (2 columns) of 1-based vertex indices

## Examples

``` r
r <- cdt_triangulate(c(0, 1, 1, 0, 0.5), c(0, 0, 1, 1, 0.5))
e <- cdt_edges(r)
plot(r$P, asp = 1)
segments(r$P[e[, 1], 1], r$P[e[, 1], 2], r$P[e[, 2], 1], r$P[e[, 2], 2])
```
