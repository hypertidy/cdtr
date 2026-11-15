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
