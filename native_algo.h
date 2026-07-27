#pragma once
// Native (Backend B) definitions of the vertex-erasure API declared in
// vertex_erasure.h. Included only under PCHS_BACKEND_NATIVE.
//
// The key simplification vs CGAL: instead of recording an ear-clipping path and
// replaying it via split_facet, we compute the convex one-ring triangulation as
// a plain triangle list (over the ring order) and splice it in one call with
// nat::Mesh::retriangulate_star.
//
// The one-ring copy is a small nat::Mesh built as a fan in ring_vertices() order.
// That order is exactly the order retriangulate_star fills the hole in, and it
// makes the fan outward-oriented (matching the dual), so primal_change sees the
// correct orientation and the spliced faces come out outward too. We convexify
// it by reusing the shared flip loop (flip_until_all_interior_edges_are_convex,
// now backend-agnostic via mesh::flip_edge).
#include "CostFunction.hpp"
#include "convex_triangulation.h"   // edge_side, flip_until_all_interior_edges_are_convex
#include "primal_change.hpp"

#include <array>
#include <utility>
#include <vector>

// Collect all neighbor vertex handles of v in one-ring order.
template <class Vertex_handle>
std::vector<Vertex_handle>
collect_neighbors(const Vertex_handle & v)
{
  std::vector<Vertex_handle> neighbors;
  auto h = v->halfedge();
  do
  {
    neighbors.push_back(h->opposite()->vertex());
    h = h->next()->opposite();
  } while(h != v->halfedge());
  return neighbors;
}

// Build the convex one-ring triangulation of v as a nat::Mesh fan (in
// ring_vertices order) + flip to convex; return the triangle list over ring
// indices and the one-ring copy (for the cost computation).
template <class Polyhedron>
std::vector<std::array<int, 3>>
build_convex_one_ring(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy)
{
  auto ring = const_cast<Polyhedron &>(dual).ring_vertices(v);
  const int k = (int)ring.size();

  std::vector<typename Polyhedron::Point_3> pts(k);
  for(int i = 0; i < k; ++i) pts[i] = ring[i]->point();

  std::vector<std::array<int, 3>> fan;
  for(int i = 1; i + 1 < k; ++i) fan.push_back({0, i, i + 1});
  one_ring_copy.build(pts, fan);

  flip_until_all_interior_edges_are_convex(one_ring_copy, 1000);

  std::vector<std::array<int, 3>> tris;
  for(auto f = one_ring_copy.facets_begin(); f != one_ring_copy.facets_end(); ++f)
  {
    auto h = f->halfedge();
    tris.push_back({h->vertex().index(),
                    h->next()->vertex().index(),
                    h->next()->next()->vertex().index()});
  }
  return tris;
}

// Cost of erasing dual vertex v + the triangulation to splice back in.
template <class Polyhedron>
std::pair<typename Polyhedron::Traits::FT, Record>
measure_vertex_erasure(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  const int /*max_degree_for_flips*/,
  const CostFunction cost_function)
{
  Polyhedron one_ring_copy;
  Record record;
  record.tris = build_convex_one_ring(dual, v, one_ring_copy);

  const auto p0 = v->point();
  auto [cost, contains_origin] = primal_change(one_ring_copy, p0, cost_function);
  (void)contains_origin;
  return {cost, record};
}

// Apply: splice the recorded triangulation into the hole left by erasing v.
template <class Polyhedron>
void erase_vertex_and_clip_ears(
  Polyhedron & dual,
  typename Polyhedron::Vertex_handle & v,
  const Record & record)
{
  dual.retriangulate_star(v, record.tris);
}
