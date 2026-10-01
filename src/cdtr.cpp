// cdtr: thin R binding to artem-ogre/CDT (MPL-2.0, header-only)
// Aimed at replacing RTriangle in hypertidy/anglr.
#include <cpp11.hpp>
#include <cpp11/matrix.hpp>
#include <limits>
#include <vector>
#include "CDT/CDT.h"

using namespace cpp11;
namespace writable = cpp11::writable;

namespace {

writable::integers unrefined_counts(const CDT::Unrefined& u) {
  writable::integers out({(int)u.shortEdgeTriangles, (int)u.circumcenterOutside,
                          (int)u.circumcenterOnVertex, (int)u.sharpFixedCorner,
                          (int)u.shortEdges, (int)u.splitVertexInvalid});
  out.names() = {"shortEdgeTriangles", "circumcenterOutside", "circumcenterOnVertex",
                 "sharpFixedCorner", "shortEdges", "splitVertexInvalid"};
  return out;
}

} // namespace

[[cpp11::register]]
list cdt_triangulate_cpp(doubles x, doubles y, integers s0, integers s1,
                         double max_area, double min_angle_deg,
                         int max_steiner, double min_edge_length,
                         bool conforming, int erase_mode, int intersect_mode,
                         bool angle_first) {
  // erase_mode: 0 = convex hull, 1 = outer (RTriangle -p semantics), 2 = outer + holes
  // intersect_mode: 0 = NotAllowed, 1 = TryResolve, 2 = DontCheck
  const R_xlen_t n = x.size();
  std::vector<CDT::V2d<double> > verts;
  verts.reserve(n);
  for (R_xlen_t i = 0; i < n; ++i) verts.push_back(CDT::V2d<double>(x[i], y[i]));

  std::vector<CDT::Edge> edges;
  edges.reserve(s0.size());
  for (R_xlen_t i = 0; i < s0.size(); ++i) {
    edges.push_back(CDT::Edge(static_cast<CDT::VertInd>(s0[i]),
                              static_cast<CDT::VertInd>(s1[i])));
  }

  // CDT requires unique vertices; keep the mapping so callers can remap
  // per-vertex attributes (PA-style) back onto the deduplicated set.
  CDT::DuplicatesInfo dup = CDT::RemoveDuplicatesAndRemapEdges(verts, edges);

  CDT::IntersectingConstraintEdges::Enum ice = CDT::IntersectingConstraintEdges::NotAllowed;
  if (intersect_mode == 1) ice = CDT::IntersectingConstraintEdges::TryResolve;
  if (intersect_mode == 2) ice = CDT::IntersectingConstraintEdges::DontCheck;

  CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, ice, 0.0);
  cdt.insertVertices(verts);
  const int n_input = static_cast<int>(verts.size());
  if (!edges.empty()) {
    if (conforming) cdt.conformToEdges(edges); else cdt.insertEdges(edges);
  }

  // Region classification before anything is erased: layer depth is the
  // number of constraint boundaries crossed from the outside (0 = outside,
  // odd = inside, even > 0 = hole), computed topologically.
  std::vector<CDT::LayerDepth> depths = cdt.calculateTriangleDepths();

  // Refinement must not spill into the super-triangle region, so hand the
  // refiner the set of triangles that will be erased.
  CDT::TriIndUSet to_erase;
  if (erase_mode == 1) to_erase = cdt.collectOuterTriangles();
  else if (erase_mode == 2) to_erase = cdt.collectOuterTrianglesAndHoles();
  else to_erase = cdt.collectSuperTriangle();

  writable::list unrefined;
  if (max_area > 0 || min_angle_deg > 0) {
    CDT::VertInd budget = max_steiner < 0 ? std::numeric_limits<CDT::VertInd>::max()
                                          : static_cast<CDT::VertInd>(max_steiner);
    // CDT refines by one criterion per call; run the requested ones in turn.
    for (int pass = 0; pass < 2; ++pass) {
      const bool do_area = angle_first ? pass == 1 : pass == 0;
      if (do_area) {
        if (max_area > 0) {
          CDT::Unrefined u = cdt.refineTriangles(budget, CDT::RefinementCriterion::LargestArea,
                                                 max_area, &to_erase, min_edge_length);
          unrefined.push_back({"area"_nm = unrefined_counts(u)});
        }
      } else if (min_angle_deg > 0) {
        CDT::Unrefined u = cdt.refineTriangles(budget, CDT::RefinementCriterion::SmallestAngle,
                                               CDT::degToRad(min_angle_deg), &to_erase, min_edge_length);
        unrefined.push_back({"angle"_nm = unrefined_counts(u)});
      }
    }
    // Steiner points change the triangle list; recompute depths.
    depths = cdt.calculateTriangleDepths();
  }

  // Depths for the triangles that survive, in surviving order (removal is
  // order-preserving in CDT::Triangulation::removeTriangles).
  std::vector<int> kept_depth;
  kept_depth.reserve(cdt.triangles.size());
  for (std::size_t i = 0; i < cdt.triangles.size(); ++i) {
    if (!to_erase.count(static_cast<CDT::TriInd>(i))) kept_depth.push_back((int)depths[i]);
  }

  cdt.finalizeTriangulation(to_erase);

  const int nv = static_cast<int>(cdt.vertices.size());
  writable::doubles_matrix<> P(nv, 2);
  for (int i = 0; i < nv; ++i) { P(i, 0) = cdt.vertices[i].x; P(i, 1) = cdt.vertices[i].y; }

  const int nt = static_cast<int>(cdt.triangles.size());
  writable::integers_matrix<> T(nt, 3);
  for (int i = 0; i < nt; ++i)
    for (int j = 0; j < 3; ++j) T(i, j) = static_cast<int>(cdt.triangles[i].vertices[j]) + 1;

  const int ns = static_cast<int>(cdt.fixedEdges.size());
  writable::integers_matrix<> S(ns, 2);
  int k = 0;
  for (CDT::EdgeUSet::const_iterator e = cdt.fixedEdges.begin(); e != cdt.fixedEdges.end(); ++e, ++k) {
    S(k, 0) = static_cast<int>(e->v1()) + 1;
    S(k, 1) = static_cast<int>(e->v2()) + 1;
  }

  writable::integers dup_map(dup.mapping.size());
  for (std::size_t i = 0; i < dup.mapping.size(); ++i) dup_map[i] = (int)dup.mapping[i] + 1;

  writable::integers depth(kept_depth.size());
  for (std::size_t i = 0; i < kept_depth.size(); ++i) depth[i] = kept_depth[i];

  return writable::list({
    "P"_nm = P,
    "T"_nm = T,
    "S"_nm = S,
    "depth"_nm = depth,
    "n_input"_nm = n_input,
    "input_map"_nm = dup_map,
    "unrefined"_nm = unrefined});
}
