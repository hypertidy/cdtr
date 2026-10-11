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
