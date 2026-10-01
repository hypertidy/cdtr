library(cdtr); library(RTriangle); library(sf); library(bench)
#S <- "/tmp/claude-0/-home-claude/d29382fe-fc08-585c-bc2f-455df9581bae/scratchpad"

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
#load(file.path(S, "anglr/data/cont_tas.rda")); load(file.path(S, "anglr/data/cad_tas.rda"))
cont_tas <- anglr::cont_tas; cad_tas <- anglr::cad_tas
cases <- list(nc = as_pslg_sf(nc),
              cont_tas = as_pslg_sf(sf::st_as_sf(cont_tas)),
              cad_tas = as_pslg_sf(sf::st_as_sf(cad_tas)))
for (nm in names(cases)) cat(nm, ": ", nrow(cases[[nm]]$P), "vertices", nrow(cases[[nm]]$S), "segments\n")

report <- function(nm, p, max_area = NULL) {
  rt <- tryCatch(if (is.null(max_area)) triangulate(p) else triangulate(p, a = max_area),
                 error = function(e) e)
  cd <- tryCatch(cdt_pslg(p, max_area = max_area), error = function(e) e)
  for (lab in c("RTriangle", "cdtr")) {
    r <- if (lab == "RTriangle") rt else cd
    if (inherits(r, "error")) { cat(sprintf("  %-9s ERROR %s\n", lab, conditionMessage(r))); next }
    ar <- tri_area(r$P, r$T); ma <- min_angle(r$P, r$T)
    cat(sprintf("  %-9s %6d verts %6d tris  area max %.3g  min angle %.2f  segs kept %s\n",
                lab, nrow(r$P), nrow(r$T), max(ar), min(ma),
                segments_preserved(p$P, p$S, r$P, r$S)))
  }
  if (!inherits(cd, "error")) cat("  cdtr depth table:", paste(names(table(cd$depth)), table(cd$depth), collapse = " "),
                                  " unrefined:", paste(unlist(cd$unrefined), collapse = ","), "\n")
  invisible(list(rt = rt, cd = cd))
}
for (nm in names(cases)) {
  p <- cases[[nm]]
  bb <- apply(p$P, 2, range); A <- prod(diff(bb))
  cat("\n==", nm, "constrained only\n"); report(nm, p)
  cat("==", nm, "max_area =", signif(A / 5000, 3), "\n"); report(nm, p, max_area = A / 5000)
}

cat("\n== timing (nc, max_area) ==\n")
p <- cases$nc; A <- prod(diff(apply(p$P, 2, range)))
print(bench::mark(RTriangle = triangulate(p, a = A / 5000),
                  cdtr = cdt_pslg(p, max_area = A / 5000), check = FALSE, iterations = 10)[, c("expression", "median", "mem_alloc")])
cat("\n== timing (cad_tas, constrained only) ==\n")
p <- cases$cad_tas
print(bench::mark(RTriangle = triangulate(p), cdtr = cdt_pslg(p), check = FALSE, iterations = 5)[, c("expression", "median", "mem_alloc")])

## crossing constraints: RTriangle inserts a crossing vertex; what does cdtr do?
cat("\n== crossing segments ==\n")
px <- pslg(P = cbind(c(0, 1, 0, 1), c(0, 1, 1, 0)), S = rbind(c(1, 2), c(3, 4)))
print(nrow(triangulate(px)$P)); print(nrow(cdt_pslg(px)$P))
## collinear overlapping segments (the shared-boundary case that crashed polymer)
cat("== collinear overlapping segments ==\n")
pc <- pslg(P = cbind(c(0, 2, 1, 3, 0, 3), c(0, 0, 0, 0, 1, 1)), S = rbind(c(1, 2), c(3, 4), c(5, 6)))
print(tryCatch(nrow(triangulate(pc)$T), error = function(e) conditionMessage(e)))
print(tryCatch(nrow(cdt_pslg(pc)$T), error = function(e) conditionMessage(e)))
