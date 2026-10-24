# cdtr (development version)

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
