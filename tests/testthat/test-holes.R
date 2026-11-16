## 4x4 square with a 2x2 hole ring, and a separate 1x1 square far away
two <- function() {
  x <- c(0, 4, 4, 0, 1, 3, 3, 1, 10, 11, 11, 10)
  y <- c(0, 0, 4, 4, 1, 1, 3, 3, 0, 0, 1, 1)
  list(x = x, y = y, s0 = c(1:3, 4, 5:7, 8, 9:11, 12),
       s1 = c(2:4, 1, 6:8, 5, 10:12, 9))
}
area <- function(r) {
  a <- r$P[r$T[, 1], ]; b <- r$P[r$T[, 2], ]; c <- r$P[r$T[, 3], ]
  sum(0.5 * abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) - (c[, 1] - a[, 1]) * (b[, 2] - a[, 2])))
}

test_that("a hole point erases the region bounded by segments", {
  d <- two()
  base <- cdt_triangulate(d$x, d$y, d$s0, d$s1)
  expect_equal(area(base), 16 + 1)
  h <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = cbind(2, 2))
  expect_equal(area(h), 16 - 4 + 1)
  ## same result as erase = "holes", with the hole given as a point
  eh <- cdt_triangulate(d$x, d$y, d$s0, d$s1, erase = "holes")
  expect_equal(area(eh), area(h))
  expect_equal(sort(h$depth), sort(base$depth[base$depth != 2L]))
})

test_that("holes can be a point in an outer ring, several points, list or data frame", {
  d <- two()
  h1 <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = cbind(0.5, 0.5))
  expect_equal(area(h1), 16 + 1 - (16 - 4))
  hm <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = rbind(c(2, 2), c(10.5, 0.5)))
  expect_equal(area(hm), 16 - 4)
  hl <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = list(x = 2, y = 2))
  hd <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = data.frame(x = 2, y = 2))
  expect_equal(hl$T, hd$T)
})

test_that("holes outside every ring, or without segments, erase nothing", {
  d <- two()
  base <- cdt_triangulate(d$x, d$y, d$s0, d$s1)
  out <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = cbind(7, 7))
  expect_equal(out$T, base$T)
  pts <- cdt_triangulate(c(0, 1, 1, 0, 0.5), c(0, 0, 1, 1, 0.5), holes = cbind(0.5, 0.2))
  expect_equal(nrow(pts$T), 4L)
})

test_that("holes work with refinement and keep the hole empty", {
  d <- two()
  r <- cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = cbind(2, 2), max_area = 0.2,
                       min_angle = 20)
  expect_equal(area(r), 16 - 4 + 1)
  expect_equal(length(cdt_locate(r$P, r$T, 2, 2)), 1L)
  expect_true(is.na(cdt_locate(r$P, r$T, 2, 2)))
  expect_true(all(r$depth != 2L))
})

test_that("cdt_pslg uses H, and explicit holes take precedence", {
  d <- two()
  p <- list(P = cbind(d$x, d$y), S = cbind(d$s0, d$s1), H = cbind(2, 2))
  expect_equal(area(cdt_pslg(p)), 16 - 4 + 1)
  expect_equal(area(cdt_pslg(p, holes = cbind(10.5, 0.5))), 16)
})

test_that("bad holes are rejected", {
  d <- two()
  expect_error(cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = cbind(NA, 1)), "finite")
  expect_error(cdt_triangulate(d$x, d$y, d$s0, d$s1, holes = matrix(1, 1, 3)), "two columns")
})
