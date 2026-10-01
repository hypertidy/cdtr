source(system.file("benchmarks", "helpers.R", package = "cdtr"))
for (nm in c("nc", "cont_tas")) {
  p <- cases[[nm]]; A <- prod(diff(apply(p$P, 2, range)))
  for (ma in c(A/5000, A/50000)) {
    rt <- triangulate(p, a = ma, q = 20)
    cd <- cdt_pslg(p, max_area = ma, min_angle = 20)
    cd2 <- cdt_pslg(p, max_area = ma, min_angle = 20, conforming = TRUE)
    f <- function(r) sprintf("%6d verts %6d tris  area max %.3g  min angle %.1f", nrow(r$P), nrow(r$T), max(tri_area(r$P, r$T)), min(min_angle(r$P, r$T)))
    cat(nm, "a =", signif(ma,3), "q = 20\n  RTriangle ", f(rt), "\n  cdtr      ", f(cd), " unrefined", paste(unlist(cd$unrefined), collapse=","), "\n  cdtr conf ", f(cd2), "\n")
    print(bench::mark(RTriangle = triangulate(p, a = ma, q = 20), cdtr = cdt_pslg(p, max_area = ma, min_angle = 20), check = FALSE, iterations = 5)[, c("expression","median")])
  }
}
