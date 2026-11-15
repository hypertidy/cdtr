// Point location in a triangle mesh: the triangle that contains each query
// point, by the shared grid of triangle buckets (grid.h) and exact
// orientation tests (the robust predicates bundled with CDT). A point on an
// edge or a vertex shared by several triangles goes to the lowest
// triangle id. Points in no triangle are NA. Triangles may be either
// orientation; degenerate (zero area) triangles never contain a point.
#include <cpp11.hpp>
#include <cpp11/matrix.hpp>
#include "CDT/CDTUtils.h"  // brings in the robust predicates
#include "grid.h"

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
integers cdt_locate_cpp(doubles_matrix<> P, integers_matrix<> T, doubles xq, doubles yq) {
  const int nv = P.nrow(), nt = T.nrow(), nq = xq.size();
  writable::integers out(nq);
  if (nq == 0) return out;
  for (int q = 0; q < nq; ++q) out[q] = NA_INTEGER;
  if (nt == 0 || nv == 0) return out;

  auto px = [&](int i) { return P(i, 0); };
  auto py = [&](int i) { return P(i, 1); };
  auto ta = [&](int t) { return T(t, 0) - 1; };
  auto tb = [&](int t) { return T(t, 1) - 1; };
  auto tc = [&](int t) { return T(t, 2) - 1; };
  cdtr::TriangleGrid g(px, py, nv, ta, tb, tc, nt);

  // orientation of each triangle, once
  std::vector<signed char> sgn(nt);
  for (int t = 0; t < nt; ++t) {
    const int a = ta(t), b = tb(t), c = tc(t);
    double o = predicates::orient2d(P(a, 0), P(a, 1), P(b, 0), P(b, 1), P(c, 0), P(c, 1));
    sgn[t] = (o > 0) - (o < 0);
  }

  for (int q = 0; q < nq; ++q) {
    const double x = xq[q], y = yq[q];
    if (std::isnan(x) || std::isnan(y)) continue;
    const int cell = g.cell_of(x, y);
    for (int k = g.start[cell]; k < g.start[cell + 1]; ++k) {
      const int t = g.tri[k];
      if (sgn[t] == 0) continue;
      if (x < g.tx0[t] || x > g.tx1[t] || y < g.ty0[t] || y > g.ty1[t]) continue;
      const int a = ta(t), b = tb(t), c = tc(t);
      const double s = sgn[t];
      if (s * predicates::orient2d(P(a, 0), P(a, 1), P(b, 0), P(b, 1), x, y) < 0) continue;
      if (s * predicates::orient2d(P(b, 0), P(b, 1), P(c, 0), P(c, 1), x, y) < 0) continue;
      if (s * predicates::orient2d(P(c, 0), P(c, 1), P(a, 0), P(a, 1), x, y) < 0) continue;
      out[q] = t + 1;
      break;
    }
  }
  return out;
}
