source(system.file("benchmarks", "helpers.R", package = "cdtr"))
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
                                  " unrefined:", fmt_unrefined(cd$unrefined), "\n")
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
