#pragma once
#include "CostFunction.hpp"
#include "geometry.h"        // pulls in the PCHS_BACKEND_* selection
#include "convex_triangulation.h"
#include "polyhedron_utils.h"

#include <array>
#include <tuple>
#include <utility>
#include <vector>

// How to re-triangulate the hole left by erasing a dual vertex. The two backends
// record it differently:
//   * CGAL: an ear-clipping path replayed via split_facet.
//   * native: the triangulation itself (triples over the ring order), spliced in
//     one call by nat::Mesh::retriangulate_star.
struct Record
{
#if defined(PCHS_BACKEND_NATIVE)
  std::vector<std::array<int, 3>> tris;
#else
  size_t start_vertex_id;
  std::vector<int> path;
#endif
};

// Collect all neighbor vertex handles of v in one-ring order.
template <class Vertex_handle>
std::vector<Vertex_handle>
collect_neighbors(const Vertex_handle & v);

// Compute the cost of erasing dual vertex v, and record how to re-triangulate
// its hole (see Record).
template <class Polyhedron>
std::pair<typename Polyhedron::Traits::FT, Record>
measure_vertex_erasure(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  const int max_degree_for_flips,
  const CostFunction cost_function = CostFunction::PRIMAL_VOLUME);

// Erase v from the dual and re-triangulate the resulting hole per `record`.
template <class Polyhedron>
void erase_vertex_and_clip_ears(
  Polyhedron & dual,
  typename Polyhedron::Vertex_handle & v,
  const Record & record);

#if !defined(PCHS_BACKEND_NATIVE)
// CGAL-only: ear-clip a one-ring copy, recording the replayable path.
template <class Polyhedron>
std::pair<size_t, std::vector<int>>
clip_ears(Polyhedron & one_ring_copy, typename Polyhedron::Halfedge_handle h0);
#endif

#if defined(PCHS_BACKEND_NATIVE)
#include "native_algo.h"     // native definitions of the above
#endif
