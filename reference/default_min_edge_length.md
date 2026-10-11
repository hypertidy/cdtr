# Default refinement edge-length floor

`min_edge_frac` times the median constraint segment length, which keeps
the refiner from cascading into sharp input corners, capped at
`area_edge_frac * sqrt(max_area)` so the floor never blocks the area
target (a triangle of area A has edges of order sqrt(A)). With no
segments only the area cap applies; with neither, 0.

## Usage

``` r
default_min_edge_length(
  x,
  y,
  s0,
  s1,
  max_area = NULL,
  min_edge_frac = 0.3,
  area_edge_frac = 0.25
)
```

## Arguments

- x, y:

  vertex coordinates

- s0, s1:

  1-based vertex indices of constraint segment start/end (may be NULL)

- max_area:

  maximum triangle area (Ruppert refinement by area), NULL for none

- min_edge_frac:

  fraction of the median constraint segment length used when
  `min_edge_length` is NULL

- area_edge_frac:

  fraction of `sqrt(max_area)` that caps the default floor

## Value

a single non-negative number, the edge-length floor that
[`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)
uses when `min_edge_length = NULL`

## Examples

``` r
x <- c(0, 1, 1, 0)
y <- c(0, 0, 1, 1)
s0 <- 1:4
s1 <- c(2:4, 1)
## 0.3 times the median segment length
default_min_edge_length(x, y, s0, s1)
#> [1] 0.3
## a small area target caps the floor at 0.25 * sqrt(max_area)
default_min_edge_length(x, y, s0, s1, max_area = 0.01)
#> [1] 0.025
## no segments and no area target: no floor
default_min_edge_length(x, y, integer(0), integer(0))
#> [1] 0
```
