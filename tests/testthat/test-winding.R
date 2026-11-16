## rings: list of list(feature, xy); returns inputs ready for cdt_triangulate()
build <- function(rings) {
  nv <- vapply(rings, function(r) nrow(r$xy), 1L)
  first <- cumsum(c(0L, nv[-length(nv)]))
  list(xy = do.call(rbind, lapply(rings, `[[`, "xy")),
       s0 = unlist(lapply(seq_along(nv), function(k) first[k] + seq_len(nv[k]))),
       s1 = unlist(lapply(seq_along(nv), function(k) first[k] + c(seq_len(nv[k])[-1], 1L))),
       owner = rep(vapply(rings, `[[`, 1, "feature"), nv))
}
inside_ring <- function(p, ring) {
  j <- nrow(ring); inside <- rep(FALSE, nrow(p))
  for (i in seq_len(nrow(ring))) {
    cross <- ((ring[i, 2] > p[, 2]) != (ring[j, 2] > p[, 2])) &
      (p[, 1] < (ring[j, 1] - ring[i, 1]) * (p[, 2] - ring[i, 2]) /
         (ring[j, 2] - ring[i, 2]) + ring[i, 1])
    inside <- xor(inside, cross); j <- i
  }
  inside
}
centres <- function(r) {
  cbind(rowMeans(matrix(r$P[t(r$T), 1], ncol = 3, byrow = TRUE)),
        rowMeans(matrix(r$P[t(r$T), 2], ncol = 3, byrow = TRUE)))
}
## coverage by even-odd point in polygon at the triangle centres
truth <- function(r, rings, nf) {
  p <- centres(r)
  sapply(seq_len(nf), function(f) {
    ins <- rep(FALSE, nrow(p))
    for (g in rings) if (g$feature == f) ins <- xor(ins, inside_ring(p, g$xy))
    ins
  })
}
sq <- function(x0, y0, x1, y1) rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1))
frame <- list(feature = 0, xy = sq(-1, -1, 7, 7))
case <- list(
  list(feature = 1, xy = sq(0, 0, 4, 4)),
  list(feature = 1, xy = sq(1, 1, 3, 3)[4:1, ]),                       # hole, clockwise
  list(feature = 2, xy = sq(2.5, 2.5, 5.5, 5.5)),
  list(feature = 3, xy = rbind(c(1, -0.5), c(3.5, 5), c(-0.5, 3))),
  list(feature = 4, xy = sq(-0.8, 0.5, 0, 3.5)),                        # shares part of an edge
  frame)

run <- function(rings, ...) {
  b <- build(rings)
  r <- cdt_triangulate(b$xy[, 1], b$xy[, 2], b$s0, b$s1, edge_sources = TRUE, ...)
  list(r = r, w = cdt_winding(r, b$owner, outside = c(-0.9, -0.9)))
}

test_that("coverage from the walk matches point in polygon", {
  o <- run(case)
  expect_false(anyNA(o$w))
  expect_equal(dim(o$w), c(nrow(o$r$T), 4L))
  expect_true(all((o$w != 0) == truth(o$r, case, 4)))
  expect_true(all(o$w %in% 0:1))
})

test_that("it still matches with refinement, and conforming mode", {
  for (args in list(list(max_area = 0.3), list(max_area = 0.5, min_angle = 25),
                    list(max_area = 0.5, conforming = TRUE))) {
    o <- do.call(run, c(list(case), args))
    expect_gt(sum(o$r$origin == "steiner"), 0)
    expect_false(anyNA(o$w))
    expect_true(all((o$w != 0) == truth(o$r, case, 4)))
  }
})

test_that("refinement does not change the coverage areas", {
  area <- function(r) {
    a <- r$P[r$T[, 1], ]; b <- r$P[r$T[, 2], ]; c <- r$P[r$T[, 3], ]
    0.5 * abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) - (c[, 1] - a[, 1]) * (b[, 2] - a[, 2]))
  }
  by_feature <- function(o) sapply(1:4, function(f) sum(area(o$r)[o$w[, f] != 0]))
  expect_equal(by_feature(run(case)), by_feature(run(case, max_area = 0.2, min_angle = 20)))
  ## feature 1 is the square less its hole, feature 2 a 3 x 3 square
  expect_equal(by_feature(run(case))[1:2], c(16 - 4, 9))
})

test_that("a clockwise outer ring has winding -1, and doubled input counts twice", {
  cw <- list(list(feature = 1, xy = sq(0, 0, 2, 2)[4:1, ]), frame)
  o <- run(cw)
  expect_true(all(o$w[o$w != 0] == -1L))
  twice <- list(list(feature = 1, xy = sq(0, 0, 2, 2)),
                list(feature = 1, xy = sq(0, 0, 2, 2)), frame)
  o2 <- run(twice)
  expect_equal(max(o2$w), 2L)
  expect_equal(sort(unique(o2$w)), c(0L, 2L))
})

test_that("random overlapping triangles agree with point in polygon", {
  set.seed(8317)
  for (rep in 1:5) {
    rings <- lapply(1:6, function(i) {
      xy <- cbind(runif(3, 0, 6), runif(3, 0, 6))
      if (sum(xy[, 1] * xy[c(2, 3, 1), 2] - xy[c(2, 3, 1), 1] * xy[, 2]) < 0) xy <- xy[3:1, ]
      list(feature = i, xy = xy)
    })
    rings <- c(rings, list(frame))
    o <- run(rings, max_area = 0.5, min_angle = 20)
    expect_false(anyNA(o$w))
    expect_true(all((o$w != 0) == truth(o$r, rings, 6)))
  }
})

test_that("segments owned by 0 are ignored", {
  o <- run(case)
  expect_false(any(o$w == 0 & truth(o$r, case, 4) == 1))
})

test_that("bad input is rejected", {
  b <- build(case)
  r_no <- cdt_triangulate(b$xy[, 1], b$xy[, 2], b$s0, b$s1)
  expect_error(cdt_winding(r_no, b$owner, c(-0.9, -0.9)), "edge_sources")
  r <- cdt_triangulate(b$xy[, 1], b$xy[, 2], b$s0, b$s1, edge_sources = TRUE)
  expect_error(cdt_winding(r, b$owner[-1], c(-0.9, -0.9)), "each input segment")
  expect_error(cdt_winding(r, b$owner, c(20, 20)), "not in the mesh")
  expect_error(cdt_winding(r, rep(0, length(b$owner)), c(-0.9, -0.9)), "at least 1")
})

test_that("triangles not connected to the start are NA with a warning", {
  ## two squares joined only by a hole-erased gap: erase = "holes" cuts the mesh
  b <- build(list(list(feature = 1, xy = sq(0, 0, 4, 4)),
                  list(feature = 1, xy = sq(1, 1, 3, 3)[4:1, ]),
                  list(feature = 2, xy = sq(1.5, 1.5, 2.5, 2.5)), frame))
  r <- cdt_triangulate(b$xy[, 1], b$xy[, 2], b$s0, b$s1, edge_sources = TRUE,
                       erase = "holes")
  expect_warning(w <- cdt_winding(r, b$owner, c(-0.9, -0.9)), "not connected")
  expect_true(anyNA(w))
})
