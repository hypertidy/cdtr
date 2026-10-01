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
#' @return list with P (vertices), T (1-based triangle indices), S (fixed edges),
#'   depth (per-triangle constraint layer depth: 0 outside, odd inside, even hole),
#'   n_input (vertices after dedupe, before Steiner), input_map (input vertex ->
#'   deduped vertex index), min_edge_length (the value used), and unrefined, a
#'   data frame with one row per refinement pass counting what CDT could not
#'   refine (short edges, sharp fixed corners, and so on)
#' @export
cdt_triangulate <- function(x, y, s0 = NULL, s1 = NULL,
                            max_area = NULL, min_angle = NULL,
                            max_steiner = Inf,
                            min_edge_length = NULL, min_edge_frac = 0.3, area_edge_frac = 0.25,
                            conforming = FALSE,
                            erase = c("outer", "hull", "holes"),
                            intersect = c("resolve", "error", "ignore"),
                            angle_first = FALSE) {
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
                             min_edge_length, conforming, erase, intersect, angle_first)
  out$min_edge_length <- min_edge_length
  out$unrefined <- unrefined_frame(out$unrefined)
  out
}

#' Default refinement edge-length floor
#'
#' `min_edge_frac` times the median constraint segment length, which keeps
#' the refiner from cascading into sharp input corners, capped at
#' `area_edge_frac * sqrt(max_area)` so the floor never blocks the area
#' target (a triangle of area A has edges of order sqrt(A)). With no
#' segments only the area cap applies; with neither, 0.
#' @inheritParams cdt_triangulate
#' @export
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
                 as.data.frame(setNames(replicate(length(cols), integer(0), simplify = FALSE), cols))))
  }
  rows <- lapply(names(u), function(nm) data.frame(criterion = nm, as.list(u[[nm]])))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Triangulate an RTriangle pslg with the CDT backend
#' @param p a pslg (list with P and S)
#' @param ... passed to cdt_triangulate
#' @export
cdt_pslg <- function(p, ...) {
  cdt_triangulate(p$P[, 1], p$P[, 2], p$S[, 1], p$S[, 2], ...)
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
#' @export
cdt_triangulate_attr <- function(x, y, s0 = NULL, s1 = NULL, PA = NULL, ...) {
  r <- cdt_triangulate(x, y, s0, s1, ...)
  if (is.null(PA)) return(r)
  PA <- as.matrix(PA)
  stopifnot(nrow(PA) == length(x))
  ## first input row for each deduplicated vertex
  first <- match(seq_len(r$n_input), r$input_map)
  base_PA <- PA[first, , drop = FALSE]
  n_out <- nrow(r$P)
  out <- matrix(NA_real_, n_out, ncol(PA), dimnames = list(NULL, colnames(PA)))
  out[seq_len(r$n_input), ] <- base_PA
  if (n_out > r$n_input) {
    base <- cdt_triangulate(x, y, s0, s1, erase = "hull", intersect = "resolve")
    idx <- seq(r$n_input + 1L, n_out)
    ## the base mesh may itself have crossing-resolution vertices beyond n_input
    bPA <- rbind(base_PA, if (nrow(base$P) > r$n_input)
      cdt_interpolate(base$P[seq_len(r$n_input), , drop = FALSE],
                      base$T, base_PA, base$P[-seq_len(r$n_input), 1], base$P[-seq_len(r$n_input), 2]))
    out[idx, ] <- cdt_interpolate(base$P, base$T, bPA, r$P[idx, 1], r$P[idx, 2])
  }
  r$PA <- out
  r
}
