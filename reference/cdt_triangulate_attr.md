# Triangulate with per-vertex attributes carried onto new vertices

[`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)
plus a `PA` matrix of per-input-vertex attributes (z, m, t ...),
returned as `$PA` aligned with the output vertices. Input vertices keep
their values; Steiner vertices get values interpolated from the
unrefined constrained triangulation of the input.

## Usage

``` r
cdt_triangulate_attr(x, y, s0 = NULL, s1 = NULL, PA = NULL, ...)
```

## Arguments

- x, y:

  vertex coordinates

- s0, s1:

  1-based vertex indices of constraint segment start/end (may be NULL)

- PA:

  numeric matrix with one row per input vertex (or NULL)

- ...:

  passed to
  [`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)

## Value

the list returned by
[`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md),
plus `PA`, a matrix with one row per output vertex and the columns of
the input `PA` (absent when `PA = NULL`)

## Examples

``` r
x <- c(0, 1, 1, 0)
y <- c(0, 0, 1, 1)
## a sloping surface: z rises with x
r <- cdt_triangulate_attr(x, y, 1:4, c(2:4, 1),
                          PA = cbind(z = 10 * x), max_area = 0.02)
## Steiner vertices carry interpolated z
head(cbind(r$P, r$PA)[r$origin == "steiner", , drop = FALSE])
#>                   z
#> [1,] 0.50 0.50  5.0
#> [2,] 0.50 0.00  5.0
#> [3,] 0.00 0.50  0.0
#> [4,] 0.25 0.25  2.5
#> [5,] 0.50 1.00  5.0
#> [6,] 1.00 0.50 10.0
all.equal(r$PA[, "z"], 10 * r$P[, 1])
#> [1] TRUE
```
