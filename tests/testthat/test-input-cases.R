## Input cases that exercise distinct paths in the bundled CDT library:
## conforming mode, vertices on constraints, overlapping constraints,
## degenerate (cocircular) input, larger inputs, and the Steiner budget.

tri_area <- function(P, T) {
  a <- P[T[, 1], , drop = FALSE]; b <- P[T[, 2], , drop = FALSE]; c <- P[T[, 3], , drop = FALSE]
  0.5 * abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) - (c[, 1] - a[, 1]) * (b[, 2] - a[, 2]))
}
edge_key <- function(m) paste(pmin(m[, 1], m[, 2]), pmax(m[, 1], m[, 2]))
seg_len <- function(P, S) sqrt(rowSums((P[S[, 1], , drop = FALSE] - P[S[, 2], , drop = FALSE])^2))

test_that("conforming mode splits a long constraint into Delaunay edges", {
  set.seed(2)
  x <- c(0, 10, runif(60, 0, 10)); y <- c(0, 3, runif(60, 0, 3))
  ## an open line bounds no region, so keep the hull
  cf <- cdt_triangulate(x, y, 1, 2, conforming = TRUE, erase = "hull")
  ## constrained mode keeps the segment whole; conforming adds Steiner points on it
  expect_equal(nrow(cdt_triangulate(x, y, 1, 2)$S), 1L)
  expect_gt(sum(cf$origin == "steiner"), 0L)
  expect_equal(sum(cf$origin == "crossing"), 0L)
  ## the pieces are collinear with the input segment and add up to its length
  steiner <- cf$P[cf$origin == "steiner", , drop = FALSE]
  expect_equal(steiner[, 2], 0.3 * steiner[, 1], tolerance = 1e-9)
  expect_equal(sum(seg_len(cf$P, cf$S)), sqrt(10^2 + 3^2), tolerance = 1e-9)
  expect_true(all(edge_key(cf$S) %in% edge_key(rbind(cf$T[, 1:2], cf$T[, 2:3], cf$T[, c(3, 1)]))))
})

test_that("conforming mode resolves crossing constraints with Steiner points", {
  cx <- cdt_triangulate(c(0, 1, 0, 1, 0.5), c(0, 1, 1, 0, 0.2), c(1, 3), c(2, 4),
                        conforming = TRUE)
  expect_equal(sum(cx$origin == "steiner"), 1L)
  expect_equal(unname(cx$P[cx$origin == "steiner", ]), c(0.5, 0.5))
  expect_equal(nrow(cx$S), 4L)
})

test_that("a vertex lying on a constraint splits it", {
  ## vertex 5 sits on the segment from 1 to 2
  x <- c(0, 2, 2, 0, 1); y <- c(0, 0, 2, 2, 0)
  r <- cdt_triangulate(x, y, 1:4, c(2:4, 1))
  expect_equal(nrow(r$P), 5L)
  expect_setequal(edge_key(r$S), edge_key(rbind(c(1, 5), c(5, 2), c(2, 3), c(3, 4), c(4, 1))))
  expect_equal(sum(tri_area(r$P, r$T)), 4)
})

test_that("overlapping collinear constraints are merged", {
  ## the inner segment 5-6 lies along the outer edge 1-2
  x <- c(0, 3, 3, 0, 1, 2); y <- c(0, 0, 3, 3, 0, 0)
  r <- cdt_triangulate(x, y, c(1, 2, 3, 4, 5), c(2, 3, 4, 1, 6))
  expect_equal(nrow(r$P), 6L)
  expect_equal(sum(seg_len(r$P, r$S)), 12)
  expect_equal(sum(tri_area(r$P, r$T)), 9)
})

test_that("intersect = 'ignore' skips the crossing check", {
  r <- cdt_triangulate(c(0, 1, 0, 1), c(0, 1, 1, 0), c(1, 3), c(2, 4), intersect = "ignore")
  expect_equal(sum(r$origin == "crossing"), 0L)
  expect_equal(nrow(r$P), 4L)
})

test_that("cocircular grid points triangulate to the full grid area", {
  ## every unit cell has four cocircular corners: exact predicates break the ties
  g <- expand.grid(x = 0:9, y = 0:9)
  r <- cdt_triangulate(g$x, g$y)
  expect_equal(nrow(r$T), 2L * 9L * 9L)
  expect_equal(tri_area(r$P, r$T), rep(0.5, 162L))
})

test_that("larger random inputs give a valid Delaunay triangulation", {
  set.seed(1)
  n <- 5000L
  x <- runif(n); y <- runif(n)
  r <- cdt_triangulate(x, y)
  expect_equal(nrow(r$P), n)
  expect_true(all(tri_area(r$P, r$T) > 0))
  ## Euler: a triangulation of n points with h on the hull has 2n - 2 - h triangles
  h <- length(grDevices::chull(x, y))
  expect_equal(nrow(r$T), 2L * n - 2L - h)
})

test_that("max_steiner caps the number of Steiner points", {
  r <- cdt_triangulate(c(0, 1, 1, 0), c(0, 0, 1, 1), 1:4, c(2:4, 1),
                       max_area = 1e-4, max_steiner = 10)
  expect_equal(sum(r$origin == "steiner"), 10L)
  expect_gt(max(tri_area(r$P, r$T)), 1e-4)
})

test_that("zero-length segments give no edge-length floor", {
  expect_equal(default_min_edge_length(c(0, 1), c(0, 1), 1, 1), 0)
})

test_that("a zero-column attribute matrix is carried through", {
  r <- cdt_triangulate_attr(c(0, 1, 1, 0), c(0, 0, 1, 1), 1:4, c(2:4, 1),
                            PA = matrix(0, 4, 0), max_area = 0.1)
  expect_equal(dim(r$PA), c(nrow(r$P), 0L))
})

test_that("attributes are interpolated across resolved crossings", {
  ## unit square boundary plus its two crossing diagonals, z = x + y, refined
  x <- c(0, 1, 1, 0); y <- c(0, 0, 1, 1)
  r <- cdt_triangulate_attr(x, y, c(1:4, 1, 2), c(2:4, 1, 3, 4), PA = cbind(z = x + y),
                            max_area = 0.02)
  expect_equal(sum(r$origin == "crossing"), 1L)
  expect_gt(sum(r$origin == "steiner"), 0L)
  expect_equal(r$PA[, "z"], r$P[, 1] + r$P[, 2], tolerance = 1e-9)
})
