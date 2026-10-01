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
