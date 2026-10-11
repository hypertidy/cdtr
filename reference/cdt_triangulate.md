# Constrained Delaunay triangulation with refinement (CDT backend)

Thin binding to the artem-ogre/CDT library. Input mirrors
RTriangle::pslg: vertices as x/y and constraint segments as index pairs
into the vertices.

## Usage

``` r
cdt_triangulate(
  x,
  y,
  s0 = NULL,
  s1 = NULL,
  max_area = NULL,
  min_angle = NULL,
  max_steiner = Inf,
  min_edge_length = NULL,
  min_edge_frac = 0.3,
  area_edge_frac = 0.25,
  conforming = FALSE,
  erase = c("outer", "hull", "holes"),
  intersect = c("resolve", "error", "ignore"),
  angle_first = FALSE,
  depth_distinct = FALSE
)
```

## Arguments

- x, y:

  vertex coordinates

- s0, s1:

  1-based vertex indices of constraint segment start/end (may be NULL)

- max_area:

  maximum triangle area (Ruppert refinement by area), NULL for none

- min_angle:

  minimum triangle angle in degrees, NULL for none

- max_steiner:

  budget of Steiner points to insert (Inf for unlimited)

- min_edge_length:

  refinement does not split edges or triangles already shorter than
  this. NULL (the default) chooses a value from the input and the
  targets (see
  [`default_min_edge_length()`](https://hypertidy.github.io/cdtr/reference/default_min_edge_length.md)):
  a fraction of the median constraint segment length, capped so that it
  never blocks the `max_area` target. This stops the refiner cascading
  into sharp input corners. Use 0 to never give up.

- min_edge_frac:

  fraction of the median constraint segment length used when
  `min_edge_length` is NULL

- area_edge_frac:

  fraction of `sqrt(max_area)` that caps the default floor

- conforming:

  if TRUE use conforming (Steiner on segments) instead of constrained

- erase:

  one of "outer" (drop triangles outside the segment-bounded region, as
  RTriangle's `-p` does), "hull" (keep the convex hull), "holes" (also
  drop even-depth regions). With no segments there is no boundary, so
  "hull" is used.

- intersect:

  how crossing constraints are handled: "resolve" inserts the crossing
  point, "error" fails, "ignore" skips the check

- angle_first:

  when both criteria are given, refine by angle before area (default is
  area then angle; the two orderings differ little)

- depth_distinct:

  also return `depth_distinct`, a per-triangle depth that counts a
  doubled (overlapping) boundary once rather than once per copy. The
  default `depth` counts per copy; the two differ only where the only
  path from outside a triangle crosses a shared, doubled boundary.

## Value

list with P (vertices), T (1-based triangle indices, wound
counter-clockwise / positive signed area), S (fixed edges), depth
(per-triangle constraint layer depth: 0 outside, odd inside, even hole;
a doubled shared boundary is counted once per copy), origin (a factor
giving each vertex's provenance: "input", "crossing" or "steiner"),
n_input (vertices after dedupe, before Steiner), n_crossing (vertices
added resolving crossing constraints), input_map (input vertex -\>
deduped vertex index), min_edge_length (the value used), and unrefined,
a data frame with one row per refinement pass counting what CDT could
not refine (short edges, sharp fixed corners, and so on). When
`depth_distinct = TRUE` the list also has depth_distinct (a doubled
boundary counted once).
