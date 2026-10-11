# Triangulate an RTriangle pslg with the CDT backend

Convenience for input in the shape of
[`RTriangle::pslg()`](https://rdrr.io/pkg/RTriangle/man/pslg.html): a
list with `P`, a two-column matrix of vertex coordinates, and `S`, a
two-column matrix of 1-based segment vertex indices.

## Usage

``` r
cdt_pslg(p, ...)
```

## Arguments

- p:

  a pslg (list with P and S)

- ...:

  passed to
  [`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)

## Value

the list returned by
[`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)

## Examples

``` r
p <- list(P = cbind(c(0, 2, 2, 1, 0), c(0, 0, 2, 1, 2)),
          S = cbind(1:5, c(2:5, 1)))
r <- cdt_pslg(p, max_area = 0.05)
nrow(r$T)
#> [1] 102
plot(r$P, asp = 1, pch = ".")
polygon(r$P[t(cbind(r$T, NA)), ])
```
