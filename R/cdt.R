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
#' @param min_edge_length refinement gives up on edges shorter than this
#' @param conforming if TRUE use conforming (Steiner on segments) instead of constrained
#' @param erase one of "hull" (keep convex hull, like RTriangle), "outer", "holes"
#' @param intersect how crossing constraints are handled: "resolve" inserts the
#'   crossing point, "error" fails, "ignore" skips the check
#' @return list with P (vertices), T (1-based triangle indices), S (fixed edges),
#'   depth (per-triangle constraint layer depth: 0 outside, odd inside, even hole),
#'   n_input (vertices after dedupe, before Steiner), input_map (input vertex ->
#'   deduped vertex index), unrefined (counts of refinements CDT could not do)
#' @export
cdt_triangulate <- function(x, y, s0 = NULL, s1 = NULL,
                            max_area = NULL, min_angle = NULL,
                            max_steiner = Inf, min_edge_length = 1e-6,
                            conforming = FALSE,
                            erase = c("outer", "hull", "holes"),
                            intersect = c("resolve", "error", "ignore"), angle_first = FALSE) {
  erase <- match(match.arg(erase), c("hull", "outer", "holes")) - 1L
  intersect <- match(match.arg(intersect), c("error", "resolve", "ignore")) - 1L
  if (is.null(s0)) s0 <- integer(0)
  if (is.null(s1)) s1 <- integer(0)
  stopifnot(length(s0) == length(s1), length(x) == length(y))
  cdt_triangulate_cpp(as.double(x), as.double(y),
                      as.integer(s0) - 1L, as.integer(s1) - 1L,
                      if (is.null(max_area)) -1 else max_area,
                      if (is.null(min_angle)) -1 else min_angle,
                      if (is.infinite(max_steiner)) -1L else as.integer(max_steiner),
                      min_edge_length, conforming, erase, intersect, angle_first)
}

#' Convert an RTriangle pslg to cdt_triangulate arguments
#' @param p a pslg
#' @param ... passed to cdt_triangulate
#' @export
cdt_pslg <- function(p, ...) {
  cdt_triangulate(p$P[, 1], p$P[, 2], p$S[, 1], p$S[, 2], ...)
}
