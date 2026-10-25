#' Constrained Delaunay triangulation with refinement (CDT backend)
#'
#' Thin binding to the artem-ogre/CDT library. Input mirrors RTriangle::pslg:
#' vertices as x/y and constraint segments as index pairs into the vertices.
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
#'   to never give up.
#' @param min_edge_frac fraction of the median constraint segment length used
#'   when `min_edge_length` is NULL
#' @param area_edge_frac fraction of `sqrt(max_area)` that caps the default floor
#' @param conforming if TRUE use conforming (Steiner on segments) instead of constrained
#' @param erase one of "outer" (drop triangles outside the segment-bounded region,
#'   as RTriangle's `-p` does), "hull" (keep the convex hull), "holes" (also drop
#'   even-depth regions). With no segments there is no boundary, so "hull" is used.
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
                            edge_sources = FALSE) {
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
                             depth_distinct, edge_sources)
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
#' @param p a pslg (list with P and S)
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
  cdt_triangulate(p$P[, 1], p$P[, 2], p$S[, 1], p$S[, 2], ...)
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
