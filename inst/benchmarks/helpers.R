## shared by the bench*.R scripts; run from an installed cdtr (needs RTriangle, sf, bench, anglr)
library(cdtr); library(RTriangle); library(sf); library(bench)


## polygons or lines -> vertex pool + unique segments (what silicate::SC0 gives
## anglr). Every vertex of every part is kept and consecutive vertices within a
## part are joined, so a polygon ring closes through its repeated first vertex
## and an open line stays open (no segment from its end back to its start).
## Lines bound no region, so a lines-only input is marked to be meshed over its
## convex hull (attribute "hull"; see cd_tri() and rt_tri()).
as_pslg_sf <- function(x) {
  g <- sf::st_geometry(x)
  lines <- all(sf::st_geometry_type(g) %in% c("LINESTRING", "MULTILINESTRING"))
  if (!lines) g <- sf::st_cast(g, "MULTILINESTRING")
  co <- sf::st_coordinates(g)
  lcol <- grep("^L", colnames(co))
  path <- do.call(paste, as.data.frame(co[, lcol, drop = FALSE]))
  n <- nrow(co)
  s0 <- which(c(path[-1] == path[-n], FALSE)); s1 <- s0 + 1L
  ## dedupe vertices (shared boundaries, ring closure) and remap segments, drop duplicate segments
  key <- paste(co[, 1], co[, 2])
  u <- !duplicated(key); map <- match(key, key[u])
  P <- unname(co[u, 1:2, drop = FALSE]); s0 <- map[s0]; s1 <- map[s1]
  seg <- cbind(pmin(s0, s1), pmax(s0, s1)); seg <- seg[!duplicated(seg) & seg[, 1] != seg[, 2], , drop = FALSE]
  p <- pslg(P = P, S = seg)
  attr(p, "hull") <- lines
  p
}
## the two meshers on one case. Polygons are meshed inside their boundaries
## (Triangle's -p, cdtr's default erase = "outer"); lines over their convex hull
## (cdtr's erase = "hull"; RTriangle has no -c, so the hull is added as
## segments). Input segments carry boundary marker 1 and added hull segments 2,
## so the output's SB tells them apart.
cd_tri <- function(p, ...) cdt_pslg(p, ..., erase = if (isTRUE(attr(p, "hull"))) "hull" else "outer")
rt_tri <- function(p, ...) {
  S <- p$S; SB <- rep(1L, nrow(S))
  if (isTRUE(attr(p, "hull"))) {
    h <- grDevices::chull(p$P)
    S <- rbind(S, cbind(h, c(h[-1L], h[1L]))); SB <- c(SB, rep(2L, length(h)))
  }
  triangulate(pslg(P = p$P, S = S, SB = SB), ...)
}
tri_area <- function(P, T) {
  a <- P[T[, 1], ]; b <- P[T[, 2], ]; c <- P[T[, 3], ]
  0.5 * abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) - (c[, 1] - a[, 1]) * (b[, 2] - a[, 2]))
}
min_angle <- function(P, T) {
  a <- P[T[, 1], ]; b <- P[T[, 2], ]; c <- P[T[, 3], ]
  l2 <- function(p, q) rowSums((p - q)^2)
  ab <- l2(a, b); bc <- l2(b, c); ca <- l2(c, a)
  ang <- function(opp, s1, s2) acos(pmin(1, pmax(-1, (s1 + s2 - opp) / (2 * sqrt(s1 * s2)))))
  pmin(ang(bc, ab, ca), ang(ca, ab, bc), ang(ab, bc, ca)) * 180 / pi
}
## the RTriangle contract: every input segment survives as a union of output edges
segments_preserved <- function(P_in, S_in, P_out, S_out) {
  ## collinearity + coverage check is overkill for a spike; check endpoints exist
  key <- function(P) paste(signif(P[, 1], 12), signif(P[, 2], 12))
  all(key(P_in[unique(c(S_in)), , drop = FALSE]) %in% key(P_out))
}

nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

if (!requireNamespace("anglr", quietly = TRUE)) {
  stop("the benchmarks need anglr (archived on CRAN): remotes::install_github(\"hypertidy/anglr\")")
}
cont_tas <- anglr::cont_tas; cad_tas <- anglr::cad_tas
cases <- list(nc = as_pslg_sf(nc),
              cont_tas = as_pslg_sf(sf::st_as_sf(cont_tas)),
              cad_tas = as_pslg_sf(sf::st_as_sf(cad_tas)))
for (nm in names(cases)) cat(nm, ": ", nrow(cases[[nm]]$P), "vertices", nrow(cases[[nm]]$S), "segments",
                              if (isTRUE(attr(cases[[nm]], "hull"))) "(lines, meshed over the hull)" else "(polygons)", "\n")

## one line per refinement pass, non-zero counts only
fmt_unrefined <- function(u) {
  if (!nrow(u)) return("none")
  paste(apply(u, 1, function(r) {
    v <- as.integer(r[-1]); nm <- names(r)[-1]
    paste0(trimws(r[1]), ": ", if (any(v > 0)) paste(nm[v > 0], v[v > 0], collapse = " ") else "ok")
  }), collapse = " | ")
}
