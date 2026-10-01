source(system.file("benchmarks", "helpers.R", package = "cdtr"))
f <- function(r) sprintf("%6d verts %6d tris  area max %.3g  min angle %.1f", nrow(r$P), nrow(r$T), max(tri_area(r$P, r$T)), min(min_angle(r$P, r$T)))
for (nm in c("nc","cont_tas")) { p <- cases[[nm]]; A <- prod(diff(apply(p$P, 2, range))); ma <- A/5000
 cat(nm, "a =", signif(ma,3), "\n")
 cat("  RTriangle a      ", f(triangulate(p, a = ma)), "\n")
 cat("  RTriangle a q20  ", f(triangulate(p, a = ma, q = 20)), "\n")
 cat("  cdtr area only   ", f(cdt_pslg(p, max_area = ma)), "\n")
 cat("  cdtr angle only  ", f(cdt_pslg(p, min_angle = 20)), "\n")
 cat("  cdtr area->angle ", f(cdt_pslg(p, max_area = ma, min_angle = 20)), "\n")
 cat("  cdtr angle->area ", f(cdt_pslg(p, max_area = ma, min_angle = 20, angle_first = TRUE)), "\n")
 print(system.time(cdt_pslg(p, max_area = ma, min_angle = 20, angle_first = TRUE))[3])
}
