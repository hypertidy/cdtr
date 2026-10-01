// Barycentric interpolation of per-vertex attributes from a triangle mesh
// onto query points, with a uniform grid over triangle bounding boxes so
// lookup is not O(n_query * n_triangle). Points that fall in no triangle
// (outside the mesh, or lost to rounding) take the nearest vertex's value.
#include <cpp11.hpp>
#include <cpp11/matrix.hpp>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
doubles_matrix<> cdt_interpolate_cpp(doubles_matrix<> P, integers_matrix<> T,
                                     doubles_matrix<> A, doubles xq, doubles yq) {
  const int nv = P.nrow(), nt = T.nrow(), k = A.ncol(), nq = xq.size();
  writable::doubles_matrix<> out(nq, k);
  if (nq == 0) return out;

  double xmin = std::numeric_limits<double>::infinity(), ymin = xmin, xmax = -xmin, ymax = -xmin;
  for (int i = 0; i < nv; ++i) {
    xmin = std::min(xmin, P(i, 0)); xmax = std::max(xmax, P(i, 0));
    ymin = std::min(ymin, P(i, 1)); ymax = std::max(ymax, P(i, 1));
  }
  const double dx = xmax - xmin, dy = ymax - ymin;
  const int ng = std::max(1, (int)std::sqrt((double)std::max(nt, 1)));
  const double cw = dx > 0 ? dx / ng : 1.0, ch = dy > 0 ? dy / ng : 1.0;
  auto cell = [&](double v, double lo, double w) {
    int c = (int)std::floor((v - lo) / w);
    return std::min(std::max(c, 0), ng - 1);
  };

  // bucket triangle indices by the grid cells their bbox covers
  std::vector<std::vector<int> > buckets(ng * ng);
  std::vector<double> tx0(nt), tx1(nt), ty0(nt), ty1(nt);
  for (int t = 0; t < nt; ++t) {
    const int a = T(t, 0) - 1, b = T(t, 1) - 1, c = T(t, 2) - 1;
    tx0[t] = std::min(P(a, 0), std::min(P(b, 0), P(c, 0)));
    tx1[t] = std::max(P(a, 0), std::max(P(b, 0), P(c, 0)));
    ty0[t] = std::min(P(a, 1), std::min(P(b, 1), P(c, 1)));
    ty1[t] = std::max(P(a, 1), std::max(P(b, 1), P(c, 1)));
    for (int i = cell(tx0[t], xmin, cw); i <= cell(tx1[t], xmin, cw); ++i)
      for (int j = cell(ty0[t], ymin, ch); j <= cell(ty1[t], ymin, ch); ++j)
        buckets[i * ng + j].push_back(t);
  }

  const double eps = 1e-9;
  for (int q = 0; q < nq; ++q) {
    const double x = xq[q], y = yq[q];
    bool found = false;
    const std::vector<int>& cand = buckets[cell(x, xmin, cw) * ng + cell(y, ymin, ch)];
    for (std::size_t ci = 0; ci < cand.size() && !found; ++ci) {
      const int t = cand[ci];
      if (x < tx0[t] - eps || x > tx1[t] + eps || y < ty0[t] - eps || y > ty1[t] + eps) continue;
      const int a = T(t, 0) - 1, b = T(t, 1) - 1, c = T(t, 2) - 1;
      const double x1 = P(a, 0), y1 = P(a, 1), x2 = P(b, 0), y2 = P(b, 1), x3 = P(c, 0), y3 = P(c, 1);
      const double det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
      if (det == 0) continue;
      const double l1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
      const double l2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
      const double l3 = 1.0 - l1 - l2;
      if (l1 < -eps || l2 < -eps || l3 < -eps) continue;
      for (int j = 0; j < k; ++j) out(q, j) = l1 * A(a, j) + l2 * A(b, j) + l3 * A(c, j);
      found = true;
    }
    if (!found) {
      int best = 0; double bd = std::numeric_limits<double>::infinity();
      for (int i = 0; i < nv; ++i) {
        const double d = (P(i, 0) - x) * (P(i, 0) - x) + (P(i, 1) - y) * (P(i, 1) - y);
        if (d < bd) { bd = d; best = i; }
      }
      for (int j = 0; j < k; ++j) out(q, j) = A(best, j);
    }
  }
  return out;
}
