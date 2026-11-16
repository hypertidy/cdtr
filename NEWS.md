# cdtr (development version)

This development version collects what will be released as cdtr 0.2.0.

## New functions

* `cdt_locate()` returns the triangle containing each query point (`NA`
  outside), using the grid of triangle buckets from `cdt_interpolate()` and
  exact orientation tests; ties on shared edges go to the lowest triangle.
  With a triangle-to-feature table this is point-in-polygon for many points.
  With `barycentric = TRUE` it also returns the barycentric coordinates of each
  point in its triangle (`b1`, `b2`, `b3`, in the order of `T`), the weights that
  interpolate any per-vertex value. The grid is one shared, flat structure
  (`src/grid.h`) used by `cdt_locate()` and `cdt_interpolate()`.

* `cdt_winding()` finds the winding number of each input feature in every
  triangle by walking the mesh with `edge_sources`, with no point-in-polygon
  step. It gives the coverage of each triangle by overlapping inputs, which
  `depth` cannot, and is unaffected by refinement.

* `cdt_quality()` returns per-triangle area, minimum and maximum angle,
  shortest and longest edge, and radius-edge ratio, to check a mesh or that
  refinement met its targets.

## New arguments and results

* `cdt_triangulate(edge_sources = TRUE)` returns `edge_sources`, one row per
  fixed edge and covering input segment (`edge`, `segment`, `forward`),
  tracked through crossings, splits and refinement. This is the
  constraint-to-input map of the mesh spec, and lets callers label regions by
  which inputs cover them, which `depth` cannot do. An input segment lying on a
  longer collinear input is no longer lost from CDT's provenance when
  refinement splits the shared edge.

* `cdt_triangulate(holes = )` takes hole points (a two-column matrix, or a list
  or data frame with `x` and `y`) and erases the region around each that is
  bounded by constraint segments, as the `H` of an `RTriangle::pslg()` does.
  `cdt_pslg()` uses `p$H` when present. A hole point in no bounded region
  erases nothing.

## Documentation

* Four vignettes: "Getting started with cdtr", "Polygons, triangle pools and
  wkpool", "A terrain surface from contour lines" and "Overlapping inputs and
  edge_sources".

* The polygon vignette shows how wkpool 0.5.0, which now uses cdtr (in
  Suggests) for the constrained triangulation behind `tri_overlay()` and
  `tri_overlay_relate()`, takes the output of `cdt_pslg()` through
  `wkpool::new_wkpool_tri()`. Use `erase = "holes"` so that hole triangles are
  dropped. polymer builds on both packages.

* The terrain vignette uses a small data file in `inst/extdata` made from the
  Tasmania 2 m DEM, The LIST, (c) State of Tasmania, CC BY 3.0 AU (see
  `inst/COPYRIGHTS`). It relates the quality of the triangles to the error and
  appearance of the surface.

* The `erase` argument documents that `"outer"` keeps the triangles inside holes.

* The documentation and README explain that refinement, by area as well as by
  angle, splits encroached constraint segments (Ruppert's rule), which gives
  more vertices and better-shaped triangles than area-only refinement in
  Triangle, and how `min_edge_length` controls it.

* The DESCRIPTION reference now uses the CRAN form `Author (year) <doi:...>`.

# cdtr 0.1.0

* First release.
* `cdt_triangulate()` wraps the 'CDT' C++ library for constrained and
  conforming Delaunay triangulation, with Ruppert refinement by maximum area
  or minimum angle, resolution of crossing constraints, and per-triangle
  topological region depth.
* Triangles are returned wound counter-clockwise (positive signed area).
* `origin` records each vertex's provenance (`input`, `crossing` or `steiner`).
* Optional `depth_distinct` counts a doubled boundary once rather than per copy.
* Helpers: `cdt_pslg()`, `cdt_triangulate_attr()`, `cdt_interpolate()`,
  `cdt_edges()` and `default_min_edge_length()`.
