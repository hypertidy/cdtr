# Polygons, triangle pools and wkpool

cdtr works on vertices and index pairs. For real polygon data, the
[wkpool](https://hypertidy.github.io/wkpool/) package reads anything
that ‘wk’ can handle, shares coincident vertices, and builds those
inputs. Its triangle pools (`wkpool_tri`) then answer the usual polygon
questions (area, centroid, DE-9IM relations, point location, dissolve)
from the triangles directly. Since wkpool 0.5.0, cdtr is used for the
constrained triangulation behind
[`tri_overlay()`](https://hypertidy.github.io/wkpool/reference/tri_overlay.html)
and
[`tri_overlay_relate()`](https://hypertidy.github.io/wkpool/reference/tri_overlay_relate.html),
and the cdtr output feeds
[`new_wkpool_tri()`](https://hypertidy.github.io/wkpool/reference/new_wkpool_tri.html)
directly. This vignette shows that route. It needs wkpool (\>= 0.5.0),
which is a suggested package.

``` r

library(cdtr)
library(wkpool)
```

## From polygons to a PSLG

Two polygons, the first with a hole, and the second overlapping its
corner.

``` r

x <- wk::wkt(c(
  "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1))",
  "POLYGON ((3.5 3.5, 6 3.5, 6 6, 3.5 6, 3.5 3.5))"
))
pool <- merge_coincident(establish_topology(x))
pslg <- as_pslg(pool)
str(pslg)
#> List of 2
#>  $ P: num [1:12, 1:2] 0 4 4 0 1 3 3 1 3.5 6 ...
#>  $ S: int [1:12, 1:2] 1 2 3 4 5 6 7 8 9 10 ...
```

[`as_pslg()`](https://hypertidy.github.io/wkpool/reference/as_pslg.html)
returns the `P` and `S` that
[`cdt_pslg()`](https://hypertidy.github.io/cdtr/reference/cdt_pslg.md)
takes.

## A triangle pool from cdtr

For a single polygon layer, triangulate and use `erase = "holes"` so the
hole triangles are dropped. (With the default `"outer"` the hole is
triangulated too, with even `depth`, and would count towards the area.)

``` r

to_tri <- function(w, ...) {
  r <- cdt_pslg(as_pslg(merge_coincident(establish_topology(w))),
                erase = "holes", ...)
  vertices <- data.frame(.vx = seq_len(nrow(r$P)), x = r$P[, 1], y = r$P[, 2])
  new_wkpool_tri(vertices, r$T[, 1], r$T[, 2], r$T[, 3])
}
a <- to_tri(x[1])
a
#> <wkpool_tri[8 triangles, 8 vertices]>
#> [1] <triangle: 6,5,2> <triangle: 8,7,4> <triangle: 7,3,4> <triangle: 5,4,1>
#> [5] <triangle: 8,4,5> <triangle: 6,2,7> <triangle: 2,3,7> <triangle: 2,5,1>
tri_measures(a)[c("area", "perimeter")]
#>   area perimeter
#> 1   12        24
```

The area is 12, the 4 by 4 square less its 2 by 2 hole.

## Questions on the mesh

Points can be located in the pool with
[`tri_locate()`](https://hypertidy.github.io/wkpool/reference/tri_locate.html),
which, like
[`cdt_locate()`](https://hypertidy.github.io/cdtr/reference/cdt_locate.md),
uses exact orientation tests.

``` r

pts <- cbind(c(0.5, 2, 3.5), c(0.5, 2, 3.5))
tri_locate(a, pts)
#> [1]  4 NA  3
```

The point in the hole, (2, 2), is in no triangle.

## Overlay between layers

[`tri_overlay()`](https://hypertidy.github.io/wkpool/reference/tri_overlay.html)
puts the boundaries of two layers through one constrained triangulation
(cdtr), so every place they meet is a shared vertex or edge. The DE-9IM
relation is then read from the overlay’s triangles, edges and vertices.

``` r

b <- to_tri(x[2])
ov <- tri_overlay(a, b)
tri_overlay_relate(ov)
#>   x y    relate
#> 1 1 1 212101212
```

The pattern `212101212` is that of two polygons whose interiors overlap
and whose boundaries cross.

## Controlling the triangulation

Arguments to
[`cdt_pslg()`](https://hypertidy.github.io/cdtr/reference/cdt_pslg.md)
such as `max_area`, `min_angle` and `conforming` go straight through, so
a pool can be refined to a quality target.

``` r

a2 <- to_tri(x[1], max_area = 0.5, min_angle = 25)
length(a2)
#> [1] 24
tri_measures(a2)$area
#> [1] 12
```

Refinement leaves the area unchanged, as it only adds Steiner vertices.
Whichever route is used, the cdtr result also has `origin` and `depth`,
which a pool built this way does not keep, so use
[`cdt_triangulate()`](https://hypertidy.github.io/cdtr/reference/cdt_triangulate.md)
directly when you need them.

For overlays of many layers, with labels for the covering features, see
the polymer package, which builds on wkpool and cdtr.
