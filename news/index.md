# Changelog

## cdtr (development version)

This development version collects what will be released as cdtr 0.2.0.

- `cdt_locate(barycentric = TRUE)` also returns the barycentric
  coordinates of each point in its triangle (`b1`, `b2`, `b3`, in the
  order of `T`), the weights that interpolate any per-vertex value.

- New
  [`cdt_winding()`](https://hypertidy.github.io/cdtr/reference/cdt_winding.md)
  finds the winding number of each input feature in every triangle by
  walking the mesh and using `edge_sources`, with no point-in-polygon
  step. It is the coverage of each triangle by overlapping inputs, which
  `depth` cannot give, and is unaffected by refinement.

- New
  [`cdt_quality()`](https://hypertidy.github.io/cdtr/reference/cdt_quality.md)
  returns per-triangle area, minimum and maximum angle, shortest and
  longest edge, and radius-edge ratio, to check a mesh or that
  refinement met its targets.

- `cdt_triangulate(holes = )` takes hole points (a two-column matrix, or
  a list or data frame with `x` and `y`) and erases the region around
  each that is bounded by constraint segments, as the `H` of an
  [`RTriangle::pslg()`](https://rdrr.io/pkg/RTriangle/man/pslg.html)
  does.
  [`cdt_pslg()`](https://hypertidy.github.io/cdtr/reference/cdt_pslg.md)
  uses `p$H` when present. This makes a drop-in replacement for
  Triangle’s holes where polygons are given as rings of segments. A hole
  point in no bounded region erases nothing.

- New vignette “Overlapping inputs and edge_sources”: triangulate
  several overlapping polygons at once and find which inputs cover each
  triangle by walking the mesh with `edge_sources`, checked against
  point-in-polygon.

- New vignette “A terrain surface from contour lines”, rebuilding a
  surface around kunanyi (Mount Wellington) from contours with
  [`cdt_triangulate_attr()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate_attr.md)
  and checking it against the elevation model. The small data file in
  `inst/extdata` is from the Tasmania 2 m DEM, The LIST, (c) State of
  Tasmania, CC BY 3.0 AU (see `inst/COPYRIGHTS`). It relates the quality
  of the triangles to the error and appearance of the surface, with a
  figure.

- Two vignettes: “Getting started with cdtr” (depth, `erase`,
  refinement, crossings, `edge_sources`, locating and interpolating) and
  “Polygons, triangle pools and wkpool”.

- wkpool 0.5.0 uses cdtr (in Suggests) for the constrained triangulation
  behind
  [`tri_overlay()`](https://hypertidy.github.io/wkpool/reference/tri_overlay.html)
  and
  [`tri_overlay_relate()`](https://hypertidy.github.io/wkpool/reference/tri_overlay_relate.html),
  and the output of
  [`cdt_pslg()`](https://hypertidy.github.io/cdtr/reference/cdt_pslg.md)
  feeds
  [`wkpool::new_wkpool_tri()`](https://hypertidy.github.io/wkpool/reference/new_wkpool_tri.html).
  The polygon vignette shows the route; use `erase = "holes"` so hole
  triangles are dropped. polymer builds on both packages.

- New
  [`cdt_locate()`](https://hypertidy.github.io/cdtr/reference/cdt_locate.md)
  returns the triangle containing each query point (`NA` outside), using
  the grid of triangle buckets from
  [`cdt_interpolate()`](https://hypertidy.github.io/cdtr/reference/cdt_interpolate.md)
  and exact orientation tests; ties on shared edges go to the lowest
  triangle. With a triangle-to-feature table this is point-in-polygon
  for many points. The grid is now one shared, flat structure
  (`src/grid.h`) used by both.

- The development version number (0.1.0.9000 and up) lets packages using
  `edge_sources` (polymer) require it.

- `cdt_triangulate(edge_sources = TRUE)` returns `edge_sources`, one row
  per fixed edge and covering input segment (`edge`, `segment`,
  `forward`), tracked through crossings, splits and refinement. This is
  the constraint-to-input map of the mesh spec, and lets callers label
  regions by which inputs cover them (winding per input), which `depth`
  cannot do. An input segment lying on a longer collinear input is no
  longer lost from CDT’s provenance when refinement splits the shared
  edge.

- DESCRIPTION reference now uses the CRAN form
  `Author (year) <doi:...>`.

## cdtr 0.1.0

CRAN release: 2026-10-21

- First release.
- [`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)
  wraps the ‘CDT’ C++ library for constrained and conforming Delaunay
  triangulation, with Ruppert refinement by maximum area or minimum
  angle, resolution of crossing constraints, and per-triangle
  topological region depth.
- Triangles are returned wound counter-clockwise (positive signed area).
- `origin` records each vertex’s provenance (`input`, `crossing` or
  `steiner`).
- Optional `depth_distinct` counts a doubled boundary once rather than
  per copy.
- Helpers:
  [`cdt_pslg()`](https://hypertidy.github.io/cdtr/reference/cdt_pslg.md),
  [`cdt_triangulate_attr()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate_attr.md),
  [`cdt_interpolate()`](https://hypertidy.github.io/cdtr/reference/cdt_interpolate.md),
  [`cdt_edges()`](https://hypertidy.github.io/cdtr/reference/cdt_edges.md)
  and
  [`default_min_edge_length()`](https://hypertidy.github.io/cdtr/reference/default_min_edge_length.md).
