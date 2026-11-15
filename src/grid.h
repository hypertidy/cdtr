// A uniform grid of triangle buckets over a mesh's extent: the point
// locator shared by cdt_interpolate() and cdt_locate(). Each triangle is
// listed in every cell its bounding box covers, in triangle order, so the
// candidates in a cell come lowest triangle first. Stored as one flat
// array with offsets per cell (no vector per cell).
#ifndef CDTR_GRID_H
#define CDTR_GRID_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace cdtr {

struct TriangleGrid {
  double xmin, ymin, cw, ch;
  int ng;
  std::vector<int> start;   // ng * ng + 1 offsets into tri
  std::vector<int> tri;     // 0-based triangle ids
  std::vector<double> tx0, tx1, ty0, ty1;   // triangle bounding boxes

  int cell(double v, double lo, double w) const {
    int c = (int) std::floor((v - lo) / w);
    return std::min(std::max(c, 0), ng - 1);
  }
  int cell_of(double x, double y) const {
    return cell(x, xmin, cw) * ng + cell(y, ymin, ch);
  }

  // px, py: vertex coordinates (nv); a, b, c: 0-based corners (nt)
  template <typename X, typename Y, typename A, typename B, typename C>
  TriangleGrid(const X& px, const Y& py, int nv, const A& a, const B& b,
               const C& c, int nt) {
    xmin = std::numeric_limits<double>::infinity();
    ymin = xmin;
    double xmax = -xmin, ymax = -xmin;
    for (int i = 0; i < nv; ++i) {
      xmin = std::min(xmin, (double) px(i)); xmax = std::max(xmax, (double) px(i));
      ymin = std::min(ymin, (double) py(i)); ymax = std::max(ymax, (double) py(i));
    }
    const double dx = xmax - xmin, dy = ymax - ymin;
    ng = std::max(1, (int) std::sqrt((double) std::max(nt, 1)));
    cw = dx > 0 ? dx / ng : 1.0;
    ch = dy > 0 ? dy / ng : 1.0;
    tx0.resize(nt); tx1.resize(nt); ty0.resize(nt); ty1.resize(nt);
    std::vector<int> i0(nt), i1(nt), j0(nt), j1(nt);
    start.assign((size_t) ng * ng + 1, 0);
    for (int t = 0; t < nt; ++t) {
      const int ia = a(t), ib = b(t), ic = c(t);
      tx0[t] = std::min((double) px(ia), std::min((double) px(ib), (double) px(ic)));
      tx1[t] = std::max((double) px(ia), std::max((double) px(ib), (double) px(ic)));
      ty0[t] = std::min((double) py(ia), std::min((double) py(ib), (double) py(ic)));
      ty1[t] = std::max((double) py(ia), std::max((double) py(ib), (double) py(ic)));
      i0[t] = cell(tx0[t], xmin, cw); i1[t] = cell(tx1[t], xmin, cw);
      j0[t] = cell(ty0[t], ymin, ch); j1[t] = cell(ty1[t], ymin, ch);
      for (int i = i0[t]; i <= i1[t]; ++i)
        for (int j = j0[t]; j <= j1[t]; ++j) start[(size_t) i * ng + j + 1]++;
    }
    for (size_t k = 1; k < start.size(); ++k) start[k] += start[k - 1];
    tri.resize(start.back());
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int t = 0; t < nt; ++t)
      for (int i = i0[t]; i <= i1[t]; ++i)
        for (int j = j0[t]; j <= j1[t]; ++j) tri[fill[(size_t) i * ng + j]++] = t;
  }
};

}  // namespace cdtr

#endif
