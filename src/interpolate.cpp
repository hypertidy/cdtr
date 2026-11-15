// Barycentric interpolation of per-vertex attributes from a triangle mesh
// onto query points, with a uniform grid over triangle bounding boxes so
// lookup is not O(n_query * n_triangle); the grid is shared with
// cdt_locate() (grid.h). Points that fall in no triangle
// (outside the mesh, or lost to rounding) take the nearest vertex's value.
#include <cpp11.hpp>
#include <cpp11/matrix.hpp>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include "grid.h"

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
doubles_matrix<> cdt_interpolate_cpp(doubles_matrix<> P, integers_matrix<> T,
                                     doubles_matrix<> A, doubles xq, doubles yq) {
  const int nv = P.nrow(), nt = T.nrow(), k = A.ncol(), nq = xq.size();
  writable::doubles_matrix<> out(nq, k);
  if (nq == 0) return out;

  auto px = [&](int i) { return P(i, 0); };
  auto py = [&](int i) { return P(i, 1); };
  auto ta = [&](int t) { return T(t, 0) - 1; };
  auto tb = [&](int t) { return T(t, 1) - 1; };
  auto tc = [&](int t) { return T(t, 2) - 1; };
  // bucket triangle indices by the grid cells their bbox covers
  cdtr::TriangleGrid g(px, py, nv, ta, tb, tc, nt);
  const std::vector<double>& tx0 = g.tx0;
  const std::vector<double>& tx1 = g.tx1;
  const std::vector<double>& ty0 = g.ty0;
  const std::vector<double>& ty1 = g.ty1;

  const double eps = 1e-9;
  for (int q = 0; q < nq; ++q) {
    const double x = xq[q], y = yq[q];
    bool found = false;
    const int cell = g.cell_of(x, y);
    for (int ci = g.start[cell]; ci < g.start[cell + 1] && !found; ++ci) {
      const int t = g.tri[ci];
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
