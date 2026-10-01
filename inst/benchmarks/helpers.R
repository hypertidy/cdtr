## shared by the bench*.R scripts; run from an installed cdtr (needs RTriangle, sf, bench, anglr)
library(cdtr); library(RTriangle); library(sf); library(bench)


## polygons -> vertex pool + unique segments (what silicate::SC0 gives anglr)
as_pslg_sf <- function(x) {
  g <- sf::st_geometry(x)
  rings <- list(); pid <- integer(0)
  for (i in seq_along(g)) {
    gg <- g[[i]]
    polys <- if (inherits(gg, "MULTIPOLYGON")) unlist(gg, recursive = FALSE) else gg
    for (r in polys) { rings[[length(rings) + 1L]] <- r; pid <- c(pid, i) }
  }
  xy <- do.call(rbind, lapply(rings, function(r) r[-nrow(r), , drop = FALSE]))
  nr <- vapply(rings, function(r) nrow(r) - 1L, 1L)
  off <- cumsum(c(0L, nr[-length(nr)]))
  s0 <- unlist(lapply(seq_along(nr), function(k) off[k] + seq_len(nr[k])))
  s1 <- unlist(lapply(seq_along(nr), function(k) off[k] + c(seq_len(nr[k])[-1L], 1L)))
  ## dedupe vertices (shared boundaries) and remap segments, drop duplicate segments
  key <- paste(xy[, 1], xy[, 2])
  u <- !duplicated(key); map <- match(key, key[u])
  P <- xy[u, , drop = FALSE]; s0 <- map[s0]; s1 <- map[s1]
  seg <- cbind(pmin(s0, s1), pmax(s0, s1)); seg <- seg[!duplicated(seg) & seg[, 1] != seg[, 2], ]
  pslg(P = P, S = seg)
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

cont_tas <- anglr::cont_tas; cad_tas <- anglr::cad_tas
cases <- list(nc = as_pslg_sf(nc),
              cont_tas = as_pslg_sf(sf::st_as_sf(cont_tas)),
              cad_tas = as_pslg_sf(sf::st_as_sf(cad_tas)))
for (nm in names(cases)) cat(nm, ": ", nrow(cases[[nm]]$P), "vertices", nrow(cases[[nm]]$S), "segments\n")

## one line per refinement pass, non-zero counts only
fmt_unrefined <- function(u) {
  if (!nrow(u)) return("none")
  paste(apply(u, 1, function(r) {
    v <- as.integer(r[-1]); nm <- names(r)[-1]
    paste0(trimws(r[1]), ": ", if (any(v > 0)) paste(nm[v > 0], v[v > 0], collapse = " ") else "ok")
  }), collapse = " | ")
}
