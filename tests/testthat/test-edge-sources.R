## rectangles as CCW rings: one row per rectangle (xmin, ymin, xmax, ymax)
rect_rings <- function(b) {
  x <- y <- s0 <- s1 <- numeric(0)
  for (i in seq_len(nrow(b))) {
    o <- length(x)
    x <- c(x, b[i, 1], b[i, 3], b[i, 3], b[i, 1])
    y <- c(y, b[i, 2], b[i, 2], b[i, 4], b[i, 4])
    s0 <- c(s0, o + 1:4)
    s1 <- c(s1, o + c(2:4, 1))
  }
  list(x = x, y = y, s0 = s0, s1 = s1)
}

## winding number of every triangle from edge_sources alone: walk from the
## hull inwards, and across each fixed edge add or subtract its sources
winding <- function(r) {
  T <- r$T
  nt <- nrow(T)
  nv <- nrow(r$P)
  a <- c(T[, 1], T[, 2], T[, 3])
  b <- c(T[, 2], T[, 3], T[, 1])
  tri <- rep(seq_len(nt), 3L)
  twin <- match(paste(b, a), paste(a, b))
  skey <- paste(r$S[, 1], r$S[, 2])
  es <- r$edge_sources
  ## change in winding crossing directed edge a->b from its left to its right
  delta <- function(i) {
    k <- match(paste(min(a[i], b[i]), max(a[i], b[i])), skey)
    if (is.na(k)) return(0)
    f <- es$forward[es$edge == k]
    same <- f == (a[i] == r$S[k, 1])
    -sum(ifelse(same, 1, -1))
  }
  w <- rep(NA_real_, nt)
  hull <- which(is.na(twin))
  w[tri[hull[1]]] <- -delta(hull[1])
  queue <- tri[hull[1]]
  bad <- 0L
  while (length(queue)) {
    t <- queue[1]; queue <- queue[-1]
    for (i in which(tri == t)) {
      if (is.na(twin[i])) {
        if (w[t] != -delta(i)) bad <- bad + 1L
        next
      }
      u <- tri[twin[i]]
      wu <- w[t] + delta(i)
      if (is.na(w[u])) {
        w[u] <- wu
        queue <- c(queue, u)
      } else if (w[u] != wu) {
        bad <- bad + 1L
      }
    }
  }
  attr(w, "inconsistent") <- bad
  w
}

test_that("edge_sources covers every fixed edge and input segment", {
  d <- rect_rings(rbind(c(0, 0, 2, 2), c(1, 1, 3, 3)))
  r <- cdt_triangulate(d$x, d$y, d$s0, d$s1, edge_sources = TRUE)
  es <- r$edge_sources
  expect_s3_class(es, "data.frame")
  expect_named(es, c("edge", "segment", "forward"))
  expect_setequal(es$edge, seq_len(nrow(r$S)))
  expect_setequal(es$segment, seq_along(d$s0))
  expect_false(is.unsorted(es$edge))
  ## forward agrees with the coordinates
  dv <- r$P[r$S[es$edge, 2], ] - r$P[r$S[es$edge, 1], ]
  ov <- cbind(d$x[d$s1[es$segment]] - d$x[d$s0[es$segment]],
              d$y[d$s1[es$segment]] - d$y[d$s0[es$segment]])
  expect_equal(rowSums(dv * ov) > 0, es$forward)
  ## absent unless asked for
  expect_null(cdt_triangulate(d$x, d$y, d$s0, d$s1)$edge_sources)
})

test_that("a short input on a longer collinear input keeps both, through refinement", {
  ## A = (0,0)-(1,0) lies on B = (0,0)-(2,0)
  x <- c(0, 1, 2, 1); y <- c(0, 0, 0, 1)
  s0 <- c(1, 1, 3, 4); s1 <- c(2, 3, 4, 1)
  for (ma in list(NULL, 0.01)) {
    r <- cdt_triangulate(x, y, s0, s1, max_area = ma, edge_sources = TRUE)
    es <- r$edge_sources
    on0 <- which(r$P[r$S[, 1], 2] == 0 & r$P[r$S[, 2], 2] == 0)
    xm <- (r$P[r$S[on0, 1], 1] + r$P[r$S[on0, 2], 1]) / 2
    for (i in seq_along(on0)) {
      seg <- sort(es$segment[es$edge == on0[i]])
      if (xm[i] < 1) expect_equal(seg, 1:2) else expect_equal(seg, 2L)
    }
  }
})

test_that("duplicated input segments are listed once per copy", {
  x <- c(0, 1, 1, 0); y <- c(0, 0, 1, 1)
  r <- cdt_triangulate(x, y, c(1:4, 1), c(2:4, 1, 2), edge_sources = TRUE)
  expect_equal(sort(as.vector(table(r$edge_sources$edge))), c(1L, 1L, 1L, 2L))
  ## zero-length segments have no rows
  r <- cdt_triangulate(x, y, c(1:4, 2), c(2:4, 1, 2), edge_sources = TRUE)
  expect_false(5L %in% r$edge_sources$segment)
})

test_that("edge_sources gives the overlap count where depth does not", {
  ## four overlapping bars around an uncovered gap [1, 2] x [1, 2]
  b <- rbind(c(0, 0, 3, 1), c(0, 2, 3, 3), c(0, 0, 1, 3), c(2, 0, 3, 3))
  d <- rect_rings(b)
  for (ma in list(NULL, 0.05)) {
    r <- cdt_triangulate(d$x, d$y, d$s0, d$s1, erase = "hull", max_area = ma,
                         edge_sources = TRUE)
    w <- winding(r)
    expect_equal(attr(w, "inconsistent"), 0L)
    cx <- (r$P[r$T[, 1], 1] + r$P[r$T[, 2], 1] + r$P[r$T[, 3], 1]) / 3
    cy <- (r$P[r$T[, 1], 2] + r$P[r$T[, 2], 2] + r$P[r$T[, 3], 2]) / 3
    count <- rowSums(outer(cx, b[, 1], ">") & outer(cx, b[, 3], "<") &
                     outer(cy, b[, 2], ">") & outer(cy, b[, 4], "<"))
    expect_equal(as.vector(w), count)
    gap <- cx > 1 & cx < 2 & cy > 1 & cy < 2
    expect_true(all(w[gap] == 0))
    expect_true(all(r$depth[gap] > 0))
  }
})
