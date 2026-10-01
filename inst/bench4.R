src <- readLines("inst/bench.R"); eval(parse(text = src[1:grep("^for \\(nm in names\\(cases\\)\\) cat", src)]))
f <- function(r) sprintf("%6d verts %6d tris  area max %.3g  min angle %.1f", nrow(r$P), nrow(r$T), max(tri_area(r$P, r$T)), min(min_angle(r$P, r$T)))
p <- cases$cont_tas; ma <- 145
cat("RTriangle q20 only ", f(triangulate(p, q = 20)), "\n")
el <- function(P, S) { d <- sqrt(rowSums((P[S[,1],] - P[S[,2],])^2)); summary(d) }
print(el(p$P, p$S))
for (mel in c(1e-6, 1, 3, 6, 12)) {
  r <- cdt_pslg(p, max_area = ma, min_angle = 20, min_edge_length = mel)
  cat(sprintf("cdtr a=145 q=20 min_edge_length=%-5g ", mel), f(r), " unrefined", paste(unlist(r$unrefined), collapse=","), "\n")
}
