# Winding number of each input feature in every triangle

For a triangulation of several overlapping inputs, find which inputs
cover each triangle. Triangles are walked from a starting point known to
be covered by nothing. Stepping across an ordinary edge changes nothing;
stepping across a fixed edge changes the winding number of each feature
covering it by plus or minus one, according to the direction of its
input segment (see `edge_sources` in
[`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)).
Only the mesh is used, with no point-in-polygon test.

## Usage

``` r
cdt_winding(x, owner, outside, n_features = max(owner))
```

## Arguments

- x:

  a result of
  [`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)
  made with `edge_sources = TRUE`

- owner:

  the feature of each input segment, an integer vector with one element
  for each segment passed as `s0`, `s1`. Zero means the segment is not
  counted, for a frame or any other construction line. Features are `1`
  to `n_features`.

- outside:

  a point `c(x, y)` in the mesh, covered by no feature

- n_features:

  the number of features, by default the largest in `owner`

## Value

an integer matrix with one row per triangle (the rows of `x$T`) and one
column per feature, giving its winding number there. Triangles the walk
cannot reach are `NA`, with a warning.

## Details

Outer rings should run counter-clockwise and holes clockwise, which
makes the winding number 1 inside the feature (a clockwise outer ring
gives -1). Overlapping and shared edges are handled, and refinement does
not matter, as `edge_sources` follows each input through splits.

The mesh needs a region covered by nothing to start from, for example a
frame around the inputs (give it feature 0, or a feature that you
ignore), and the walk must be able to reach every triangle. Use
`erase = "hull"`, or the default `"outer"` with a frame. With
`erase = "holes"` the triangles are not connected through the erased
regions, so some are not reached.

## Examples

``` r
## two overlapping squares, in a frame (feature 0)
sq <- function(x0, y0, x1, y1) rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1))
rings <- list(sq(0, 0, 2, 2), sq(1, 1, 3, 3), sq(-1, -1, 4, 4))
feature <- c(1, 2, 0)
n <- vapply(rings, nrow, 1L)
first <- cumsum(c(0L, n[-length(n)]))
xy <- do.call(rbind, rings)
s0 <- unlist(lapply(seq_along(n), function(k) first[k] + 1:4))
s1 <- unlist(lapply(seq_along(n), function(k) first[k] + c(2:4, 1L)))
r <- cdt_triangulate(xy[, 1], xy[, 2], s0, s1, max_area = 0.5,
                     edge_sources = TRUE)
w <- cdt_winding(r, owner = rep(feature, n), outside = c(-0.5, -0.5))
## triangles by the features covering them
table(feature_1 = w[, 1], feature_2 = w[, 2])
#>          feature_2
#> feature_1  0  1
#>         0 48  6
#>         1  6  2
```
