// cdtr: thin R binding to artem-ogre/CDT (MPL-2.0, header-only)
// Aimed at replacing RTriangle in hypertidy/anglr.
#include <cpp11.hpp>
#include <cpp11/matrix.hpp>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <utility>
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
                         bool angle_first, bool want_depth_distinct,
                         bool want_edge_sources) {
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
  // Vertex count after the input is in (includes the super-triangle); the
  // difference against the post-edge count is how many vertices edge handling
  // added, i.e. crossing-constraint resolution points (or conforming points).
  const std::size_t nv_after_insert = cdt.vertices.size();
  if (!edges.empty()) {
    if (conforming) cdt.conformToEdges(edges); else cdt.insertEdges(edges);
  }
  const int n_crossing = static_cast<int>(cdt.vertices.size() - nv_after_insert);

  // Edge provenance: which input segments cover each fixed edge. CDT records
  // the originals of a fixed edge in pieceToOriginals only when the edge is
  // not itself an input edge, so an input edge that is also covered by a
  // longer collinear input edge lists only the longer one; and if that piece
  // is later split by refinement, the shorter input is lost. Adding the edge
  // itself to its own originals before refinement keeps every input.
  std::unordered_map<CDT::Edge, std::vector<int> > input_of;
  if (want_edge_sources) {
    for (std::size_t i = 0; i < edges.size(); ++i) {
      if (edges[i].v1() != edges[i].v2()) input_of[edges[i]].push_back((int)i);
    }
    for (CDT::EdgeUSet::const_iterator e = cdt.fixedEdges.begin(); e != cdt.fixedEdges.end(); ++e) {
      std::unordered_map<CDT::Edge, CDT::EdgeVec>::iterator it = cdt.pieceToOriginals.find(*e);
      if (it == cdt.pieceToOriginals.end()) continue;
      const CDT::Edge self(CDT::VertInd(e->v1() - CDT::nSuperTriVerts),
                           CDT::VertInd(e->v2() - CDT::nSuperTriVerts));
      if (input_of.count(self) &&
          std::find(it->second.begin(), it->second.end(), *e) == it->second.end()) {
        it->second.push_back(*e);
      }
    }
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

  // Optional depth that counts a doubled (overlapping) boundary once rather
  // than once per copy: clear the overlap counts and redo the same topological
  // peel, then restore them.
  std::vector<CDT::LayerDepth> depths_distinct;
  if (want_depth_distinct) {
    auto saved_overlap = cdt.overlapCount;
    cdt.overlapCount.clear();
    depths_distinct = cdt.calculateTriangleDepths();
    cdt.overlapCount = saved_overlap;
  }

  // Depths for the triangles that survive, in surviving order (removal is
  // order-preserving in CDT::Triangulation::removeTriangles).
  std::vector<int> kept_depth, kept_depth_distinct;
  kept_depth.reserve(cdt.triangles.size());
  for (std::size_t i = 0; i < cdt.triangles.size(); ++i) {
    if (!to_erase.count(static_cast<CDT::TriInd>(i))) {
      kept_depth.push_back((int)depths[i]);
      if (want_depth_distinct) kept_depth_distinct.push_back((int)depths_distinct[i]);
    }
  }

  cdt.finalizeTriangulation(to_erase);

  const int nv = static_cast<int>(cdt.vertices.size());
  writable::doubles_matrix<> P(nv, 2);
  for (int i = 0; i < nv; ++i) { P(i, 0) = cdt.vertices[i].x; P(i, 1) = cdt.vertices[i].y; }

  const int nt = static_cast<int>(cdt.triangles.size());
  writable::integers_matrix<> T(nt, 3);
  for (int i = 0; i < nt; ++i) {
    int a = static_cast<int>(cdt.triangles[i].vertices[0]);
    int b = static_cast<int>(cdt.triangles[i].vertices[1]);
    int c = static_cast<int>(cdt.triangles[i].vertices[2]);
    // Enforce counter-clockwise winding (positive signed area).
    const double area2 =
      (cdt.vertices[b].x - cdt.vertices[a].x) * (cdt.vertices[c].y - cdt.vertices[a].y) -
      (cdt.vertices[c].x - cdt.vertices[a].x) * (cdt.vertices[b].y - cdt.vertices[a].y);
    if (area2 < 0) std::swap(b, c);
    T(i, 0) = a + 1; T(i, 1) = b + 1; T(i, 2) = c + 1;
  }

  const int ns = static_cast<int>(cdt.fixedEdges.size());
  writable::integers_matrix<> S(ns, 2);
  int k = 0;
  for (CDT::EdgeUSet::const_iterator e = cdt.fixedEdges.begin(); e != cdt.fixedEdges.end(); ++e, ++k) {
    S(k, 0) = static_cast<int>(e->v1()) + 1;
    S(k, 1) = static_cast<int>(e->v2()) + 1;
  }

  // One row per (fixed edge, covering input segment). `forward` is TRUE when
  // the input segment runs the same way as the edge as listed in S.
  writable::integers es_edge, es_segment;
  writable::logicals es_forward;
  if (want_edge_sources) {
    k = 0;
    for (CDT::EdgeUSet::const_iterator e = cdt.fixedEdges.begin(); e != cdt.fixedEdges.end(); ++e, ++k) {
      CDT::EdgeVec originals;
      std::unordered_map<CDT::Edge, CDT::EdgeVec>::const_iterator it = cdt.pieceToOriginals.find(*e);
      if (it != cdt.pieceToOriginals.end()) originals = it->second;
      if (std::find(originals.begin(), originals.end(), *e) == originals.end() && input_of.count(*e)) {
        originals.push_back(*e);
      }
      const double dx = cdt.vertices[e->v2()].x - cdt.vertices[e->v1()].x;
      const double dy = cdt.vertices[e->v2()].y - cdt.vertices[e->v1()].y;
      for (std::size_t j = 0; j < originals.size(); ++j) {
        std::unordered_map<CDT::Edge, std::vector<int> >::const_iterator in = input_of.find(originals[j]);
        if (in == input_of.end()) continue;
        for (std::size_t m = 0; m < in->second.size(); ++m) {
          const int seg = in->second[m];
          const CDT::VertInd a = dup.mapping[s0[seg]];
          const CDT::VertInd b = dup.mapping[s1[seg]];
          const double ox = cdt.vertices[b].x - cdt.vertices[a].x;
          const double oy = cdt.vertices[b].y - cdt.vertices[a].y;
          es_edge.push_back(k + 1);
          es_segment.push_back(seg + 1);
          es_forward.push_back(dx * ox + dy * oy > 0 ? TRUE : FALSE);
        }
      }
    }
  }

  writable::integers dup_map(dup.mapping.size());
  for (std::size_t i = 0; i < dup.mapping.size(); ++i) dup_map[i] = (int)dup.mapping[i] + 1;

  writable::integers depth(kept_depth.size());
  for (std::size_t i = 0; i < kept_depth.size(); ++i) depth[i] = kept_depth[i];

  writable::list out({
    "P"_nm = P,
    "T"_nm = T,
    "S"_nm = S,
    "depth"_nm = depth,
    "n_input"_nm = n_input,
    "n_crossing"_nm = n_crossing,
    "input_map"_nm = dup_map,
    "unrefined"_nm = unrefined});

  if (want_depth_distinct) {
    writable::integers depth_distinct(kept_depth_distinct.size());
    for (std::size_t i = 0; i < kept_depth_distinct.size(); ++i)
      depth_distinct[i] = kept_depth_distinct[i];
    out.push_back({"depth_distinct"_nm = depth_distinct});
  }

  if (want_edge_sources) {
    writable::list es({"edge"_nm = es_edge, "segment"_nm = es_segment,
                       "forward"_nm = es_forward});
    out.push_back({"edge_sources"_nm = es});
  }

  return out;
}
