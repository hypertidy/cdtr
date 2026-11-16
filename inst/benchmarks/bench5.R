## The numbers behind the README's "Refinement protects the constraints": with
## an area bound of bounding box area / 5000, how many Steiner vertices cdtr and
## Triangle add, how many of cdtr's lie on the input constraints, how many
## triangles have a minimum angle under 5 degrees, and what the min_edge_length
## floor does to cdtr's vertex count and to the area target.
source(system.file("benchmarks", "helpers.R", package = "cdtr"))
under5 <- function(r) mean(min_angle(r$P, r$T) < 5)
met <- function(r, ma) mean(tri_area(r$P, r$T) <= ma * (1 + 1e-9))
for (nm in names(cases)) {
  p <- cases[[nm]]; n0 <- nrow(p$P); ma <- prod(diff(apply(p$P, 2, range))) / 5000
  cat(sprintf("\n== %s: max_area %.4g (sqrt %.3g)\n", nm, ma, sqrt(ma)))
  ta <- rt_tri(p, a = ma); tq <- rt_tri(p, a = ma, q = 20)
  for (lab in c("Triangle a", "Triangle a q=20")) {
    r <- if (lab == "Triangle a") ta else tq
    st <- seq_len(nrow(r$P))[-seq_len(n0)]
    on <- mean(st %in% unique(c(r$S[r$SB == 1L, ])))
    cat(sprintf("  %-26s added %6d  on constraints %5.1f%%  under 5 deg %5.1f%%  area met %6.2f%%\n",
                lab, length(st), 100 * on, 100 * under5(r), 100 * met(r, ma)))
  }
  for (mel in list(NULL, 0, sqrt(ma) / 2, sqrt(ma))) {
    r <- cd_tri(p, max_area = ma, min_edge_length = mel)
    st <- which(r$origin == "steiner")
    on <- mean(st %in% unique(c(r$S)))
    cat(sprintf("  cdtr floor %-15s added %6d  on constraints %5.1f%%  under 5 deg %5.1f%%  area met %6.2f%%  unrefined: %s\n",
                format(signif(r$min_edge_length, 3)), length(st), 100 * on, 100 * under5(r), 100 * met(r, ma),
                fmt_unrefined(r$unrefined)))
  }
}
