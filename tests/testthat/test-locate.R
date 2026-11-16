test_that("cdt_locate() finds the containing triangle", {
  P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
  T <- rbind(c(1, 2, 3), c(1, 3, 4))
  expect_identical(
    cdt_locate(P, T, c(0.75, 0.25, 0.5, 2, NA), c(0.25, 0.75, 0.5, 2, 0)),
    c(1L, 2L, 1L, NA, NA)
  )
  ## vertices and outer edges are inside (closed triangles)
  expect_identical(cdt_locate(P, T, c(0, 1, 0.5), c(1, 0, 0)), c(2L, 1L, 1L))
  ## clockwise triangles work too
  expect_identical(cdt_locate(P, T[, 3:1], 0.75, 0.25), 1L)
  expect_identical(cdt_locate(P, T, numeric(0), numeric(0)), integer(0))
})

test_that("cdt_locate() agrees with a brute force test on a refined mesh", {
  set.seed(1)
  ang <- seq(0, 2 * pi, length.out = 41)[-41]
  x <- cos(ang); y <- sin(ang)
  n <- length(x)
  tr <- cdt_triangulate(x, y, s0 = seq_len(n), s1 = c(seq_len(n)[-1], 1L), max_area = 0.001)
  P <- tr$P; T <- tr$T
  xq <- runif(2000, -1.1, 1.1); yq <- runif(2000, -1.1, 1.1)
  got <- cdt_locate(P, T, xq, yq)
  inside <- function(t, px, py) {
    a <- P[T[t, 1], ]; b <- P[T[t, 2], ]; c <- P[T[t, 3], ]
    o <- function(u, v) (v[1] - u[1]) * (py - u[2]) - (v[2] - u[2]) * (px - u[1])
    o(a, b) >= 0 & o(b, c) >= 0 & o(c, a) >= 0
  }
  brute <- vapply(seq_along(xq), function(i) {
    hit <- which(vapply(seq_len(nrow(T)), inside, logical(1), px = xq[i], py = yq[i]))
    if (length(hit)) hit[1] else NA_integer_
  }, 1L)
  expect_identical(got, brute)
})

test_that("cdt_locate() checks its triangle indices", {
  P <- cbind(c(0, 1, 0), c(0, 0, 1))
  expect_error(cdt_locate(P, rbind(c(1, 2, 4)), 0, 0), "index rows")
})

test_that("barycentric coordinates reproduce the point and interpolate", {
  set.seed(2954)
  xy <- cbind(runif(40), runif(40))
  r <- cdt_triangulate(xy[, 1], xy[, 2])
  q <- cbind(runif(200, -0.2, 1.2), runif(200, -0.2, 1.2))
  loc <- cdt_locate(r$P, r$T, q[, 1], q[, 2], barycentric = TRUE)
  expect_equal(names(loc), c("triangle", "b1", "b2", "b3"))
  expect_equal(loc$triangle, cdt_locate(r$P, r$T, q[, 1], q[, 2]))
  inside <- !is.na(loc$triangle)
  expect_true(any(inside) && any(!inside))
  expect_true(all(is.na(loc[!inside, c("b1", "b2", "b3")])))
  b <- as.matrix(loc[inside, c("b1", "b2", "b3")])
  expect_true(all(b >= 0))
  expect_equal(unname(rowSums(b)), rep(1, nrow(b)))
  tr <- r$T[loc$triangle[inside], ]
  back <- cbind(rowSums(b * matrix(r$P[tr, 1], ncol = 3)),
                rowSums(b * matrix(r$P[tr, 2], ncol = 3)))
  expect_equal(back, q[inside, ], ignore_attr = TRUE)
  ## same as the interpolation of a per-vertex value
  z <- sin(3 * r$P[, 1]) + r$P[, 2]
  zb <- rowSums(b * matrix(z[tr], ncol = 3))
  zi <- cdt_interpolate(r$P, r$T, z, q[inside, 1], q[inside, 2])[, 1]
  expect_equal(unname(zb), zi, tolerance = 1e-12)
})

test_that("vertices and edges give exact weights, and either winding works", {
  P <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
  T <- rbind(c(1, 2, 3), c(1, 3, 4))
  v <- cdt_locate(P, T, c(0, 1, 1), c(0, 0, 1), barycentric = TRUE)
  expect_equal(unname(as.matrix(v[c("b1", "b2", "b3")])),
               rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)))
  e <- cdt_locate(P, T, 0.5, 0, barycentric = TRUE)   # midpoint of edge 1-2
  expect_equal(unlist(e[c("b1", "b2", "b3")]), c(b1 = 0.5, b2 = 0.5, b3 = 0))
  cw <- cdt_locate(P, T[, c(1, 3, 2)], 0.75, 0.25, barycentric = TRUE)
  ccw <- cdt_locate(P, T, 0.75, 0.25, barycentric = TRUE)
  expect_equal(unlist(cw[c("b1", "b3", "b2")]), unlist(ccw[c("b1", "b2", "b3")]),
               ignore_attr = TRUE)
  expect_equal(nrow(cdt_locate(P, T, numeric(0), numeric(0), barycentric = TRUE)), 0L)
  expect_true(is.na(cdt_locate(P, T, NA_real_, 0.5, barycentric = TRUE)$b1))
})
