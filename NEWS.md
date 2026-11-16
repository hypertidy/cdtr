# cdtr (development version)

This development version collects what will be released as cdtr 0.2.0.

* New `cdt_quality()` returns per-triangle area, minimum and maximum angle,
  shortest and longest edge, and radius-edge ratio, to check a mesh or that
  refinement met its targets.

* `cdt_triangulate(holes = )` takes hole points (a two-column matrix, or a
  list or data frame with `x` and `y`) and erases the region around each that
  is bounded by constraint segments, as the `H` of an `RTriangle::pslg()` does.
  `cdt_pslg()` uses `p$H` when present. This makes a drop-in replacement for
  Triangle's holes where polygons are given as rings of segments. A hole point
  in no bounded region erases nothing.

* New vignette "A terrain surface from contour lines", rebuilding a surface
  around kunanyi (Mount Wellington) from contours with `cdt_triangulate_attr()`
  and checking it against the elevation model. The small data file in
  `inst/extdata` is from the Tasmania 2 m DEM, The LIST, (c) State of Tasmania,
  CC BY 3.0 AU (see `inst/COPYRIGHTS`).

* Two vignettes: "Getting started with cdtr" (depth, `erase`, refinement,
  crossings, `edge_sources`, locating and interpolating) and "Polygons,
  triangle pools and wkpool".

* wkpool 0.5.0 uses cdtr (in Suggests) for the constrained triangulation
  behind `tri_overlay()` and `tri_overlay_relate()`, and the output of
  `cdt_pslg()` feeds `wkpool::new_wkpool_tri()`. The polygon vignette shows
  the route; use `erase = "holes"` so hole triangles are dropped. polymer
  builds on both packages.

* New `cdt_locate()` returns the triangle containing each query point (`NA`
  outside), using the grid of triangle buckets from `cdt_interpolate()` and
  exact orientation tests; ties on shared edges go to the lowest triangle.
  With a triangle-to-feature table this is point-in-polygon for many points.
  The grid is now one shared, flat structure (`src/grid.h`) used by both.

* The development version number (0.1.0.9000 and up) lets packages using
  `edge_sources` (polymer) require it.

* `cdt_triangulate(edge_sources = TRUE)` returns `edge_sources`, one row per
  fixed edge and covering input segment (`edge`, `segment`, `forward`),
  tracked through crossings, splits and refinement. This is the
  constraint-to-input map of the mesh spec, and lets callers label regions
  by which inputs cover them (winding per input), which `depth` cannot do.
  An input segment lying on a longer collinear input is no longer lost from
  CDT's provenance when refinement splits the shared edge.
* DESCRIPTION reference now uses the CRAN form `Author (year) <doi:...>`.

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
