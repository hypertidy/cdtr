## unit square with a square hole: 8 vertices, two closed rings
sq <- function() {
  x <- c(0, 1, 1, 0, 0.25, 0.75, 0.75, 0.25)
  y <- c(0, 0, 1, 1, 0.25, 0.25, 0.75, 0.75)
  s0 <- c(1, 2, 3, 4, 5, 6, 7, 8); s1 <- c(2, 3, 4, 1, 6, 7, 8, 5)
  list(x = x, y = y, s0 = s0, s1 = s1)
}
tri_area <- function(P, T) {
  a <- P[T[, 1], , drop = FALSE]; b <- P[T[, 2], , drop = FALSE]; c <- P[T[, 3], , drop = FALSE]
  0.5 * abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) - (c[, 1] - a[, 1]) * (b[, 2] - a[, 2]))
}

test_that("constrained triangulation keeps segments and classifies depth", {
  d <- sq()
  r <- cdt_triangulate(d$x, d$y, d$s0, d$s1)
  expect_equal(nrow(r$P), 8L)
  expect_equal(nrow(r$S), 8L)
  ## outer erase keeps the hole triangles (depth 2) and the ring (depth 1)
  expect_setequal(unique(r$depth), c(1L, 2L))
  expect_equal(sum(tri_area(r$P, r$T)), 1)
  expect_equal(sum(tri_area(r$P, r$T)[r$depth == 1L]), 1 - 0.25)
  ## all fixed edges appear as triangle edges
  tri_edges <- rbind(r$T[, 1:2], r$T[, 2:3], r$T[, c(3, 1)])
  key <- function(m) paste(pmin(m[, 1], m[, 2]), pmax(m[, 1], m[, 2]))
  expect_true(all(key(r$S) %in% key(tri_edges)))
})

test_that("erase modes differ as expected", {
  d <- sq()
  holes <- cdt_triangulate(d$x, d$y, d$s0, d$s1, erase = "holes")
  expect_true(all(holes$depth == 1L))
  expect_equal(sum(tri_area(holes$P, holes$T)), 0.75)
  hull <- cdt_triangulate(d$x, d$y, d$s0, d$s1, erase = "hull")
  expect_equal(sum(tri_area(hull$P, hull$T)), 1)
})

test_that("area refinement honours the bound and reports min_edge_length", {
  d <- sq()
  r <- cdt_triangulate(d$x, d$y, d$s0, d$s1, max_area = 0.01)
  expect_true(all(tri_area(r$P, r$T) <= 0.01 + 1e-12))
  expect_gt(nrow(r$P), 8L)
  expect_equal(r$n_input, 8L)
  expect_s3_class(r$unrefined, "data.frame")
  expect_equal(r$unrefined$criterion, "area")
  ## default floor: min(0.3 * median segment length (0.75), 0.25 * sqrt(max_area))
  expect_equal(r$min_edge_length, 0.25 * sqrt(0.01))
  expect_equal(cdt_triangulate(d$x, d$y, d$s0, d$s1, min_angle = 20)$min_edge_length, 0.3 * 0.75)
  expect_equal(cdt_triangulate(d$x, d$y, max_area = 0.01)$min_edge_length, 0.25 * sqrt(0.01))
  r0 <- cdt_triangulate(d$x, d$y, d$s0, d$s1, max_area = 0.01, min_edge_length = 0)
  expect_equal(r0$min_edge_length, 0)
  both <- cdt_triangulate(d$x, d$y, d$s0, d$s1, max_area = 0.01, min_angle = 20)
  expect_equal(both$unrefined$criterion, c("area", "angle"))
  expect_equal(cdt_triangulate(d$x, d$y, d$s0, d$s1, max_area = 0.01, min_angle = 20,
                               angle_first = TRUE)$unrefined$criterion, c("angle", "area"))
})

test_that("no refinement gives an empty unrefined frame and zero floor", {
  d <- sq()
  r <- cdt_triangulate(d$x, d$y, d$s0, d$s1)
  expect_equal(nrow(r$unrefined), 0L)
  expect_equal(r$min_edge_length, 0)
})

test_that("duplicate input vertices are removed and mapped", {
  r <- cdt_triangulate(c(0, 1, 0, 1, 0), c(0, 0, 1, 1, 0))
  expect_equal(nrow(r$P), 4L)
  expect_equal(r$input_map, c(1L, 2L, 3L, 4L, 1L))
})

test_that("crossing constraints are resolved by inserting the crossing", {
  r <- cdt_triangulate(c(0, 1, 0, 1), c(0, 1, 1, 0), c(1, 3), c(2, 4))
  expect_equal(nrow(r$P), 5L)
  expect_error(cdt_triangulate(c(0, 1, 0, 1), c(0, 1, 1, 0), c(1, 3), c(2, 4), intersect = "error"))
})

test_that("cdt_pslg accepts an RTriangle-style list", {
  d <- sq()
  p <- list(P = cbind(d$x, d$y), S = cbind(d$s0, d$s1))
  expect_equal(cdt_pslg(p)$T, cdt_triangulate(d$x, d$y, d$s0, d$s1)$T)
})
