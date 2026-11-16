#' Constrained Delaunay triangulation with refinement (CDT backend)
#'
#' Thin binding to the artem-ogre/CDT library. Input mirrors RTriangle::pslg:
#' vertices as x/y and constraint segments as index pairs into the vertices.
#'
#' @details
#' Refinement follows Ruppert's algorithm, which protects the constraint
#' segments: a segment is split whenever a vertex lies inside the circle on its
#' diameter (it is "encroached"), and this applies to refinement by `max_area`
#' alone as well as by `min_angle`. A mesh refined by `max_area` therefore has
#' more vertices along the constraints, and better-shaped triangles, than one
#' from a refiner that only splits large triangles, such as Triangle with `a`
#' but no `q`. With input that has many close or crossing segments, such as
#' contour lines, the difference can be a factor of two or more in the number of
#' vertices.
#'
#' `min_edge_length` is the control for this: edges shorter than it are not
#' split. Raising it gives fewer vertices. A floor of about half of
#' `sqrt(max_area)` is a reasonable first try, and a floor near `sqrt(max_area)`
#' may leave a few triangles over the target. The `unrefined` component of the
#' result counts what was left undone, so check it and the largest area after
#' raising the floor.
#'
#' @param x,y vertex coordinates
#' @param s0,s1 1-based vertex indices of constraint segment start/end (may be NULL)
#' @param max_area maximum triangle area (Ruppert refinement by area), NULL for none
#' @param min_angle minimum triangle angle in degrees, NULL for none
#' @param max_steiner budget of Steiner points to insert (Inf for unlimited)
#' @param min_edge_length refinement does not split edges or triangles already
#'   shorter than this. NULL (the default) chooses a value from the input and
#'   the targets (see [default_min_edge_length()]): a fraction of the median
#'   constraint segment length, capped so that it never blocks the `max_area`
#'   target. This stops the refiner cascading into sharp input corners. Use 0
#'   to never give up. A larger value gives a coarser mesh, since it also stops
#'   the splitting of encroached constraint segments (see Details).
#' @param min_edge_frac fraction of the median constraint segment length used
#'   when `min_edge_length` is NULL
#' @param area_edge_frac fraction of `sqrt(max_area)` that caps the default floor
#' @param conforming if TRUE use conforming (Steiner on segments) instead of constrained
#' @param erase one of "outer" (drop triangles outside the segment-bounded region,
#'   as RTriangle's `-p` does), "hull" (keep the convex hull), "holes" (also drop
#'   even-depth regions). With no segments there is no boundary, so "hull" is used.
#'   Note that `"outer"` keeps the triangles inside holes (they have even
#'   `depth`), so use `"holes"` when the triangles are to represent polygon
#'   area, for example to build a triangle pool or sum areas.
#' @param holes optional hole points, a two-column matrix or a list or data
#'   frame with `x` and `y`. The region around each point that is bounded by
#'   constraint segments is erased, as with the `H` of an `RTriangle::pslg()`.
#'   This is in addition to `erase`. A point in no bounded region (outside every
#'   ring, or with no segments) erases nothing, and a point exactly on a segment
#'   is not supported. Holes are applied after the topological `depth` is
#'   computed, so `depth` still describes the unerased mesh.
#' @param intersect how crossing constraints are handled: "resolve" inserts the
#'   crossing point, "error" fails, "ignore" skips the check
#' @param angle_first when both criteria are given, refine by angle before area
#'   (default is area then angle; the two orderings differ little)
#' @param depth_distinct also return `depth_distinct`, a per-triangle depth that
#'   counts a doubled (overlapping) boundary once rather than once per copy. The
#'   default `depth` counts per copy; the two differ only where the only path
#'   from outside a triangle crosses a shared, doubled boundary.
#' @param edge_sources also return `edge_sources`, a data frame with one row
#'   per fixed edge and input segment covering it (see Value)
#' @return list with P (vertices), T (1-based triangle indices, wound
#'   counter-clockwise / positive signed area), S (fixed edges), depth
#'   (per-triangle constraint layer depth: 0 outside, odd inside, even hole;
#'   a doubled shared boundary is counted once per copy), origin (a factor
#'   giving each vertex's provenance: "input", "crossing" or "steiner"),
#'   n_input (vertices after dedupe, before Steiner), n_crossing (vertices
#'   added resolving crossing constraints), input_map (input vertex -> deduped
#'   vertex index), min_edge_length (the value used), and unrefined, a data
#'   frame with one row per refinement pass counting what CDT could not refine
#'   (short edges, sharp fixed corners, and so on). When `depth_distinct = TRUE`
#'   the list also has depth_distinct (a doubled boundary counted once).
#'   When `edge_sources = TRUE` the list also has edge_sources, a data frame
#'   with columns `edge` (row in S), `segment` (position of the input segment
#'   in `s0`, `s1`) and `forward` (TRUE when the input segment runs the same
#'   way as the S row, from its first to its second vertex). Every input
#'   segment covering an edge has its own row, through crossings, splits and
#'   refinement, so an edge shared by k input segments (duplicates included)
#'   has k rows. Zero-length input segments have no rows.
#' @export
#' @examples
#' ## unit square with a square hole: two closed rings of 4 vertices
#' x <- c(0, 1, 1, 0, 0.25, 0.75, 0.75, 0.25)
#' y <- c(0, 0, 1, 1, 0.25, 0.25, 0.75, 0.75)
#' s0 <- c(1, 2, 3, 4, 5, 6, 7, 8)
#' s1 <- c(2, 3, 4, 1, 6, 7, 8, 5)
#'
#' ## constrained triangulation, hole triangles kept with even depth
#' r <- cdt_triangulate(x, y, s0, s1)
#' table(r$depth)
#'
#' ## drop the hole and refine by area and minimum angle
#' r2 <- cdt_triangulate(x, y, s0, s1, erase = "holes",
#'                       max_area = 0.005, min_angle = 25)
#' table(r2$origin)
#' plot(r2$P, asp = 1, pch = ".")
#' polygon(r2$P[t(cbind(r2$T, NA)), ], col = "grey90")
#'
#' ## erase the region around a point instead (like RTriangle's pslg H)
#' r5 <- cdt_triangulate(x, y, s0, s1, holes = cbind(0.5, 0.5))
#' table(r5$depth)
#'
#' ## which input segments make up each fixed edge
#' r4 <- cdt_triangulate(x, y, s0, s1, max_area = 0.05, edge_sources = TRUE)
#' head(r4$edge_sources)
#'
#' ## points only: Delaunay triangulation of the convex hull
#' set.seed(1)
#' r3 <- cdt_triangulate(runif(30), runif(30))
#' nrow(r3$T)
cdt_triangulate <- function(x, y, s0 = NULL, s1 = NULL,
                            max_area = NULL, min_angle = NULL,
                            max_steiner = Inf,
                            min_edge_length = NULL, min_edge_frac = 0.3, area_edge_frac = 0.25,
                            conforming = FALSE,
                            erase = c("outer", "hull", "holes"),
                            intersect = c("resolve", "error", "ignore"),
                            angle_first = FALSE, depth_distinct = FALSE,
                            edge_sources = FALSE, holes = NULL) {
  hxy <- hole_points_xy(holes)
  if (is.null(s0)) s0 <- integer(0)
  if (is.null(s1)) s1 <- integer(0)
  erase <- match.arg(erase)
  ## with no constraints there is no boundary to erase outside of
  if (length(s0) == 0L && erase != "hull") erase <- "hull"
  erase <- match(erase, c("hull", "outer", "holes")) - 1L
  intersect <- match(match.arg(intersect), c("error", "resolve", "ignore")) - 1L
  stopifnot(length(s0) == length(s1), length(x) == length(y))
  x <- as.double(x); y <- as.double(y)
  s0 <- as.integer(s0); s1 <- as.integer(s1)
  refining <- !is.null(max_area) || !is.null(min_angle)
  if (is.null(min_edge_length)) {
    min_edge_length <- if (refining) default_min_edge_length(x, y, s0, s1, max_area, min_edge_frac, area_edge_frac) else 0
  }
  out <- cdt_triangulate_cpp(x, y, s0 - 1L, s1 - 1L,
                             if (is.null(max_area)) -1 else max_area,
                             if (is.null(min_angle)) -1 else min_angle,
                             if (is.infinite(max_steiner)) -1L else as.integer(max_steiner),
                             min_edge_length, conforming, erase, intersect, angle_first,
                             depth_distinct, edge_sources, hxy$x, hxy$y)
  if (edge_sources) {
    es <- out$edge_sources
    out$edge_sources <- data.frame(edge = es$edge, segment = es$segment,
                                   forward = es$forward)
  }
  out$min_edge_length <- min_edge_length
  out$unrefined <- unrefined_frame(out$unrefined)
  out$origin <- vertex_origin(nrow(out$P), out$n_input, out$n_crossing, conforming)
  out
}

hole_points_xy <- function(holes) {
  if (is.null(holes)) return(list(x = double(0), y = double(0)))
  if (is.data.frame(holes) || (is.list(holes) && !is.null(holes$x))) {
    xy <- list(x = holes$x, y = holes$y)
  } else {
    holes <- as.matrix(holes)
    if (ncol(holes) != 2L) stop("`holes` must have two columns (x, y)", call. = FALSE)
    xy <- list(x = holes[, 1], y = holes[, 2])
  }
  xy <- lapply(xy, as.double)
  if (anyNA(unlist(xy)) || !all(is.finite(unlist(xy)))) {
    stop("`holes` must be finite", call. = FALSE)
  }
  xy
}

## Per-vertex provenance (meshspec `origin`): input vertices first, then
## vertices added while handling constraints, then refinement Steiner points.
## Constraint-handling points are crossing-resolution points in the usual case;
## in conforming mode they are on-segment subdivisions, so count as Steiner.
vertex_origin <- function(nv, n_input, n_crossing, conforming) {
  n_steiner <- nv - n_input - n_crossing
  mid <- if (conforming) "steiner" else "crossing"
  factor(rep(c("input", mid, "steiner"), c(n_input, n_crossing, n_steiner)),
         levels = c("input", "crossing", "steiner"))
}

#' Default refinement edge-length floor
#'
#' `min_edge_frac` times the median constraint segment length, which keeps
#' the refiner from cascading into sharp input corners, capped at
#' `area_edge_frac * sqrt(max_area)` so the floor never blocks the area
#' target (a triangle of area A has edges of order sqrt(A)). With no
#' segments only the area cap applies; with neither, 0.
#' @inheritParams cdt_triangulate
#' @return a single non-negative number, the edge-length floor that
#'   [cdt_triangulate()] uses when `min_edge_length = NULL`
#' @export
#' @examples
#' x <- c(0, 1, 1, 0)
#' y <- c(0, 0, 1, 1)
#' s0 <- 1:4
#' s1 <- c(2:4, 1)
#' ## 0.3 times the median segment length
#' default_min_edge_length(x, y, s0, s1)
#' ## a small area target caps the floor at 0.25 * sqrt(max_area)
#' default_min_edge_length(x, y, s0, s1, max_area = 0.01)
#' ## no segments and no area target: no floor
#' default_min_edge_length(x, y, integer(0), integer(0))
default_min_edge_length <- function(x, y, s0, s1, max_area = NULL,
                                    min_edge_frac = 0.3, area_edge_frac = 0.25) {
  seg <- Inf
  if (length(s0) > 0L) {
    len <- sqrt((x[s0] - x[s1])^2 + (y[s0] - y[s1])^2)
    len <- len[len > 0]
    if (length(len) > 0L) seg <- min_edge_frac * stats::median(len)
  }
  cap <- if (is.null(max_area)) Inf else area_edge_frac * sqrt(max_area)
  out <- min(seg, cap)
  if (is.infinite(out)) 0 else out
}

unrefined_frame <- function(u) {
  cols <- c("shortEdgeTriangles", "circumcenterOutside", "circumcenterOnVertex",
            "sharpFixedCorner", "shortEdges", "splitVertexInvalid")
  if (length(u) == 0L) {
    return(cbind(data.frame(criterion = character(0)),
                 as.data.frame(stats::setNames(replicate(length(cols), integer(0), simplify = FALSE), cols))))
  }
  rows <- lapply(names(u), function(nm) data.frame(criterion = nm, as.list(u[[nm]])))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Triangulate an RTriangle pslg with the CDT backend
#'
#' Convenience for input in the shape of `RTriangle::pslg()`: a list with `P`,
#' a two-column matrix of vertex coordinates, and `S`, a two-column matrix of
#' 1-based segment vertex indices.
#'
#' @param p a pslg (list with P and S, and optionally H, a two-column matrix
#'   of hole points, passed as `holes` unless `holes` is given in `...`)
#' @param ... passed to [cdt_triangulate()]
#' @return the list returned by [cdt_triangulate()]
#' @export
#' @examples
#' p <- list(P = cbind(c(0, 2, 2, 1, 0), c(0, 0, 2, 1, 2)),
#'           S = cbind(1:5, c(2:5, 1)))
#' r <- cdt_pslg(p, max_area = 0.05)
#' nrow(r$T)
#' plot(r$P, asp = 1, pch = ".")
#' polygon(r$P[t(cbind(r$T, NA)), ])
cdt_pslg <- function(p, ...) {
  args <- list(...)
  if (!is.null(p$H) && length(p$H) > 0L && is.null(args$holes)) args$holes <- p$H
  do.call(cdt_triangulate, c(list(p$P[, 1], p$P[, 2], p$S[, 1], p$S[, 2]), args))
}

#' Unique edges of a triangulation
#'
#' The index pairs of every triangle edge once, for plotting with
#' [graphics::segments()] or for building an edge table.
#'
#' @param x result of [cdt_triangulate()]
#' @return integer matrix (2 columns) of 1-based vertex indices
#' @export
#' @examples
#' r <- cdt_triangulate(c(0, 1, 1, 0, 0.5), c(0, 0, 1, 1, 0.5))
#' e <- cdt_edges(r)
#' plot(r$P, asp = 1)
#' segments(r$P[e[, 1], 1], r$P[e[, 1], 2], r$P[e[, 2], 1], r$P[e[, 2], 2])
cdt_edges <- function(x) {
  T <- x$T
  e <- rbind(T[, 1:2], T[, 2:3], T[, c(3, 1)])
  e <- cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]))
  unique(e)
}

#' Interpolate per-vertex attributes onto points
#'
#' Linear (barycentric) interpolation of a matrix of per-vertex values from a
#' triangle mesh onto query points. Points in no triangle take the nearest
#' vertex's value. Used to carry attributes such as z onto Steiner vertices
#' after refinement, matching the `PA` behaviour of RTriangle.
#'
#' @param P vertex coordinates (n x 2)
#' @param T triangle vertex indices (m x 3, 1-based)
#' @param A per-vertex attribute matrix (n x k) or vector
#' @param xq,yq query coordinates
#' @return a length(xq) x k matrix (column names kept from `A`)
#' @export
#' @examples
#' ## two triangles on the unit square, z = x + y at the corners
#' P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
#' T <- rbind(c(1, 2, 3), c(1, 3, 4))
#' A <- cbind(z = P[, 1] + P[, 2])
#' ## a linear field is reproduced exactly inside; (2, 2) is outside
#' cdt_interpolate(P, T, A, xq = c(0.5, 0.25, 2), yq = c(0.5, 0.1, 2))
cdt_interpolate <- function(P, T, A, xq, yq) {
  A <- as.matrix(A); storage.mode(A) <- "double"
  P <- as.matrix(P); storage.mode(P) <- "double"
  T <- as.matrix(T); storage.mode(T) <- "integer"
  stopifnot(nrow(A) == nrow(P), ncol(P) == 2L, ncol(T) == 3L, length(xq) == length(yq))
  out <- cdt_interpolate_cpp(P, T, A, as.double(xq), as.double(yq))
  colnames(out) <- colnames(A)
  out
}

#' Locate points in a triangle mesh
#'
#' Find the triangle that contains each query point. This is the point
#' locator inside [cdt_interpolate()] (a uniform grid of triangle buckets),
#' with exact orientation tests from the bundled robust predicates, so a
#' point is never lost to rounding between two triangles. With a triangle to
#' feature table it answers point-in-polygon for many points at once.
#'
#' A point on an edge or vertex shared by several triangles is given to the
#' lowest triangle index that contains it. Points outside every triangle, and
#' `NA` points, give `NA`. Triangles may be in either orientation, and
#' degenerate (zero area) triangles contain no point.
#'
#' With `barycentric = TRUE` the barycentric coordinates of each point in its
#' triangle are returned as well. They are the weights of the triangle's three
#' vertices (in the order given in `T`) that reproduce the point, so
#' `b1 * P[T[i, 1], ] + b2 * P[T[i, 2], ] + b3 * P[T[i, 3], ]` is the point,
#' and any per-vertex value is interpolated by the same weights, as
#' [cdt_interpolate()] does. The weights are not negative and sum to 1; they are
#' computed in double precision, with rounding error clamped, so a point on an
#' edge has a weight of exactly 0 for the opposite vertex.
#'
#' @inheritParams cdt_interpolate
#' @param barycentric if `TRUE`, return a data frame with the barycentric
#'   coordinates, rather than only the triangle index
#' @return by default an integer vector, `length(xq)`, of 1-based row indices
#'   into `T` (`NA` outside the mesh). With `barycentric = TRUE`, a data frame
#'   with one row per query point and columns `triangle` (the same index),
#'   and `b1`, `b2`, `b3`, which are `NA` where `triangle` is.
#' @export
#' @examples
#' P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
#' T <- rbind(c(1, 2, 3), c(1, 3, 4))
#' ## (0.5, 0.5) is on the shared diagonal: the lower triangle id wins
#' cdt_locate(P, T, xq = c(0.75, 0.25, 0.5, 2), yq = c(0.25, 0.75, 0.5, 2))
#'
#' ## the weights of each triangle's vertices
#' loc <- cdt_locate(P, T, xq = c(0.75, 0.25, 0.5), yq = c(0.25, 0.75, 0.5),
#'                   barycentric = TRUE)
#' loc
#' ## a value known at the vertices, interpolated with those weights
#' z <- c(0, 10, 20, 10)
#' rowSums(as.matrix(loc[c("b1", "b2", "b3")]) * matrix(z[T[loc$triangle, ]], ncol = 3))
cdt_locate <- function(P, T, xq, yq, barycentric = FALSE) {
  P <- as.matrix(P); storage.mode(P) <- "double"
  T <- as.matrix(T); storage.mode(T) <- "integer"
  stopifnot(ncol(P) == 2L, ncol(T) == 3L, length(xq) == length(yq))
  if (nrow(T) > 0L && (anyNA(T) || min(T) < 1L || max(T) > nrow(P))) {
    stop("`T` must index rows of `P`", call. = FALSE)
  }
  xq <- as.double(xq); yq <- as.double(yq)
  tri <- cdt_locate_cpp(P, T, xq, yq)
  if (!barycentric) return(tri)
  na <- rep(NA_real_, length(tri))
  out <- data.frame(triangle = tri, b1 = na, b2 = na, b3 = na)
  i <- which(!is.na(tri))
  if (length(i) > 0L) {
    a <- P[T[tri[i], 1L], , drop = FALSE]
    b <- P[T[tri[i], 2L], , drop = FALSE]
    c <- P[T[tri[i], 3L], , drop = FALSE]
    cross <- function(p, q, r) (q[, 1] - p[, 1]) * (r[, 2] - p[, 2]) -
      (q[, 2] - p[, 2]) * (r[, 1] - p[, 1])
    pt <- cbind(xq[i], yq[i])
    total <- cross(a, b, c)
    w <- cbind(cross(pt, b, c), cross(a, pt, c), cross(a, b, pt)) / total
    ## located by exact tests, so any negative weight is rounding
    w[w < 0] <- 0
    w <- w / rowSums(w)
    out$b1[i] <- w[, 1]; out$b2[i] <- w[, 2]; out$b3[i] <- w[, 3]
  }
  out
}

#' Triangulate with per-vertex attributes carried onto new vertices
#'
#' `cdt_triangulate()` plus a `PA` matrix of per-input-vertex attributes
#' (z, m, t ...), returned as `$PA` aligned with the output vertices. Input
#' vertices keep their values; Steiner vertices get values interpolated from
#' the unrefined constrained triangulation of the input.
#'
#' @inheritParams cdt_triangulate
#' @param PA numeric matrix with one row per input vertex (or NULL)
#' @param ... passed to [cdt_triangulate()]
#' @return the list returned by [cdt_triangulate()], plus `PA`, a matrix with
#'   one row per output vertex and the columns of the input `PA` (absent when
#'   `PA = NULL`)
#' @export
#' @examples
#' x <- c(0, 1, 1, 0)
#' y <- c(0, 0, 1, 1)
#' ## a sloping surface: z rises with x
#' r <- cdt_triangulate_attr(x, y, 1:4, c(2:4, 1),
#'                           PA = cbind(z = 10 * x), max_area = 0.02)
#' ## Steiner vertices carry interpolated z
#' head(cbind(r$P, r$PA)[r$origin == "steiner", , drop = FALSE])
#' all.equal(r$PA[, "z"], 10 * r$P[, 1])
cdt_triangulate_attr <- function(x, y, s0 = NULL, s1 = NULL, PA = NULL, ...) {
  r <- cdt_triangulate(x, y, s0, s1, ...)
  if (is.null(PA)) return(r)
  PA <- as.matrix(PA)
  stopifnot(nrow(PA) == length(x))
  n_out <- nrow(r$P)
  if (ncol(PA) == 0L) {
    r$PA <- matrix(0, n_out, 0L)
    return(r)
  }
  ## first input row for each deduplicated vertex
  first <- match(seq_len(r$n_input), r$input_map)
  base_PA <- PA[first, , drop = FALSE]
  out <- matrix(NA_real_, n_out, ncol(PA), dimnames = list(NULL, colnames(PA)))
  out[seq_len(r$n_input), ] <- base_PA
  if (n_out > r$n_input) {
    base <- cdt_triangulate(x, y, s0, s1, erase = "hull", intersect = "resolve")
    idx <- seq(r$n_input + 1L, n_out)
    ## the base mesh may itself have crossing-resolution vertices beyond n_input
    bPA <- base_PA
    if (nrow(base$P) > r$n_input) {
      extra <- seq(r$n_input + 1L, nrow(base$P))
      bPA <- rbind(base_PA, cdt_interpolate(base$P[seq_len(r$n_input), , drop = FALSE],
                                            base$T, base_PA, base$P[extra, 1], base$P[extra, 2]))
    }
    out[idx, ] <- cdt_interpolate(base$P, base$T, bPA, r$P[idx, 1], r$P[idx, 2])
  }
  r$PA <- out
  r
}

#' Triangle quality measures
#'
#' Per-triangle area, angles, edge lengths and radius-edge ratio, for judging
#' a mesh or checking that refinement met its targets.
#'
#' @param x a result of [cdt_triangulate()] (a list with `P` and `T`), or a
#'   two-column matrix of vertex coordinates, in which case `T` is needed
#' @param T triangle vertex indices (m x 3, 1-based), when `x` is a matrix
#' @return a data frame with one row per triangle and columns `area` (always
#'   non-negative), `min_angle` and `max_angle` (degrees), `min_edge` and
#'   `max_edge`, and `radius_edge` (circumradius over shortest edge; the
#'   Ruppert bound on it is about 1.41 for a 20.7 degree minimum angle, and
#'   `Inf` for a degenerate triangle). Use [summary()] on the result for the
#'   distributions.
#' @export
#' @examples
#' ## an irregular polygon with a sharp corner, refined by area and angle
#' x <- c(0, 2, 2, 1, 0)
#' y <- c(0, 0, 2, 0.4, 2)
#' r <- cdt_triangulate(x, y, 1:5, c(2:5, 1), max_area = 0.05, min_angle = 25)
#' q <- cdt_quality(r)
#' summary(q[c("area", "min_angle", "radius_edge")])
#' ## both refinement targets were met
#' c(max_area = max(q$area), min_angle = min(q$min_angle))
cdt_quality <- function(x, T = NULL) {
  if (is.list(x) && !is.null(x$P) && !is.null(x$T)) {
    P <- x$P; T <- x$T
  } else {
    P <- x
    if (is.null(T)) stop("`T` is needed when `x` is a matrix", call. = FALSE)
  }
  P <- as.matrix(P); T <- as.matrix(T)
  stopifnot(ncol(P) == 2L, ncol(T) == 3L)
  if (nrow(T) > 0L && (anyNA(T) || min(T) < 1L || max(T) > nrow(P))) {
    stop("`T` must index rows of `P`", call. = FALSE)
  }
  ax <- P[T[, 1], 1]; ay <- P[T[, 1], 2]
  bx <- P[T[, 2], 1]; by <- P[T[, 2], 2]
  cx <- P[T[, 3], 1]; cy <- P[T[, 3], 2]
  a <- sqrt((bx - cx)^2 + (by - cy)^2)  # opposite vertex 1
  b <- sqrt((ax - cx)^2 + (ay - cy)^2)
  c <- sqrt((ax - bx)^2 + (ay - by)^2)
  area <- 0.5 * abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))
  ang <- function(opp, s1, s2) {
    cs <- (s1^2 + s2^2 - opp^2) / (2 * s1 * s2)
    acos(pmin(1, pmax(-1, cs))) * 180 / pi
  }
  angles <- cbind(ang(a, b, c), ang(b, a, c), ang(c, a, b))
  min_edge <- pmin(a, b, c)
  max_edge <- pmax(a, b, c)
  circumradius <- a * b * c / (4 * area)
  data.frame(area = area,
             min_angle = suppressWarnings(apply(angles, 1, min)),
             max_angle = suppressWarnings(apply(angles, 1, max)),
             min_edge = min_edge, max_edge = max_edge,
             radius_edge = circumradius / min_edge)
}

#' Winding number of each input feature in every triangle
#'
#' For a triangulation of several overlapping inputs, find which inputs cover
#' each triangle. Triangles are walked from a starting point known to be
#' covered by nothing. Stepping across an ordinary edge changes nothing;
#' stepping across a fixed edge changes the winding number of each feature
#' covering it by plus or minus one, according to the direction of its input
#' segment (see `edge_sources` in [cdt_triangulate()]). Only the mesh is used,
#' with no point-in-polygon test.
#'
#' Outer rings should run counter-clockwise and holes clockwise, which makes
#' the winding number 1 inside the feature (a clockwise outer ring gives -1).
#' Overlapping and shared edges are handled, and refinement does not matter, as
#' `edge_sources` follows each input through splits.
#'
#' The mesh needs a region covered by nothing to start from, for example a
#' frame around the inputs (give it feature 0, or a feature that you ignore),
#' and the walk must be able to reach every triangle. Use `erase = "hull"`, or
#' the default `"outer"` with a frame. With `erase = "holes"` the triangles are
#' not connected through the erased regions, so some are not reached.
#'
#' @param x a result of [cdt_triangulate()] made with `edge_sources = TRUE`
#' @param owner the feature of each input segment, an integer vector with one
#'   element for each segment passed as `s0`, `s1`. Zero means the segment
#'   is not counted, for a frame or any other construction line. Features are
#'   `1` to `n_features`.
#' @param outside a point `c(x, y)` in the mesh, covered by no feature
#' @param n_features the number of features, by default the largest in `owner`
#' @return an integer matrix with one row per triangle (the rows of `x$T`) and
#'   one column per feature, giving its winding number there. Triangles the walk
#'   cannot reach are `NA`, with a warning.
#' @export
#' @examples
#' ## two overlapping squares, in a frame (feature 0)
#' sq <- function(x0, y0, x1, y1) rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1))
#' rings <- list(sq(0, 0, 2, 2), sq(1, 1, 3, 3), sq(-1, -1, 4, 4))
#' feature <- c(1, 2, 0)
#' n <- vapply(rings, nrow, 1L)
#' first <- cumsum(c(0L, n[-length(n)]))
#' xy <- do.call(rbind, rings)
#' s0 <- unlist(lapply(seq_along(n), function(k) first[k] + 1:4))
#' s1 <- unlist(lapply(seq_along(n), function(k) first[k] + c(2:4, 1L)))
#' r <- cdt_triangulate(xy[, 1], xy[, 2], s0, s1, max_area = 0.5,
#'                      edge_sources = TRUE)
#' w <- cdt_winding(r, owner = rep(feature, n), outside = c(-0.5, -0.5))
#' ## triangles by the features covering them
#' table(feature_1 = w[, 1], feature_2 = w[, 2])
cdt_winding <- function(x, owner, outside, n_features = max(owner)) {
  if (is.null(x$edge_sources)) {
    stop("`x` needs `edge_sources`: call cdt_triangulate(edge_sources = TRUE)", call. = FALSE)
  }
  es <- x$edge_sources
  owner <- as.integer(owner)
  if (anyNA(owner) || any(owner < 0L)) stop("`owner` must be non-negative integers", call. = FALSE)
  if (nrow(es) > 0L && length(owner) < max(es$segment)) {
    stop("`owner` needs one element for each input segment", call. = FALSE)
  }
  if (n_features < 1L) stop("`n_features` must be at least 1", call. = FALSE)
  if (any(owner > n_features)) stop("`owner` has values above `n_features`", call. = FALSE)
  P <- x$P; T <- x$T; S <- x$S
  nt <- nrow(T); ns <- nrow(S); nv <- nrow(P)
  w <- matrix(NA_integer_, nt, n_features)
  if (nt == 0L) return(w)

  ## winding change of each feature on stepping across an S row from its right
  ## to its left; a segment running against the S row counts the other way
  delta <- matrix(0L, ns, n_features)
  f <- owner[es$segment]
  keep <- f > 0L
  if (any(keep)) {
    cell <- (f[keep] - 1L) * ns + es$edge[keep]
    sums <- rowsum(ifelse(es$forward[keep], 1L, -1L), cell)
    delta[as.integer(rownames(sums))] <- as.integer(sums)
  }

  ## triangle edge k joins vertices k and k + 1 (mod 3); a counter-clockwise
  ## triangle lies to the left of that edge as directed
  a <- as.vector(T[, c(1L, 2L, 3L)])
  b <- as.vector(T[, c(2L, 3L, 1L)])
  tri <- rep(seq_len(nt), 3L)
  k <- rep(1:3, each = nt)
  key <- pmin(a, b) * (nv + 1) + pmax(a, b)
  ord <- order(key)
  nbr_tri <- nbr_edge <- matrix(NA_integer_, nt, 3L)
  i <- which(key[ord][-length(ord)] == key[ord][-1L])
  p <- ord[i]; q <- ord[i + 1L]
  nbr_tri[cbind(tri[p], k[p])] <- tri[q]; nbr_edge[cbind(tri[p], k[p])] <- k[q]
  nbr_tri[cbind(tri[q], k[q])] <- tri[p]; nbr_edge[cbind(tri[q], k[q])] <- k[p]

  s_key <- pmin(S[, 1], S[, 2]) * (nv + 1) + pmax(S[, 1], S[, 2])
  s_row <- matrix(match(key, s_key), nt, 3L)

  start <- cdt_locate(P, T, outside[1L], outside[2L])
  if (is.na(start)) stop("`outside` is not in the mesh", call. = FALSE)
  w[start, ] <- 0L
  queue <- integer(nt); queue[1L] <- start
  head <- 1L; tail <- 1L
  while (head <= tail) {
    t <- queue[head]; head <- head + 1L
    for (e in 1:3) {
      u <- nbr_tri[t, e]
      if (is.na(u) || !anyNA(w[u, ])) next
      sr <- s_row[t, e]
      if (is.na(sr)) {
        w[u, ] <- w[t, ]
      } else {
        ## u is on the left of the S row if its own edge runs the same way
        ue <- nbr_edge[t, e]
        same_way <- T[u, ue] == S[sr, 1L]
        w[u, ] <- w[t, ] + (if (same_way) 1L else -1L) * delta[sr, ]
      }
      tail <- tail + 1L; queue[tail] <- u
    }
  }
  if (anyNA(w[, 1L])) warning("some triangles are not connected to `outside`", call. = FALSE)
  w
}
