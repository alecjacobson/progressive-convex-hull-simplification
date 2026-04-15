#pragma once
#include "CostFunction.hpp"
#include "convex_triangulation.h"
#include "polyhedron_utils.h"

#include <tuple>
#include <utility>
#include <vector>

// Encodes the ear-clipping triangulation order for a one-ring so it can be
// replayed on the actual dual without re-triangulating.
struct Record
{
  size_t start_vertex_id;
  std::vector<int> path;
};

// Ear-clip one_ring_copy (a triangulated disk with open boundary) starting
// from h0. Records the clipping order as a path (walks between ears) and the
// starting vertex id so it can be replayed with erase_vertex_and_clip_ears.
template <class Polyhedron>
std::pair<size_t, std::vector<int>>
clip_ears(Polyhedron & one_ring_copy, typename Polyhedron::Halfedge_handle h0);


// Apply a previously recorded ear-clipping path to the actual dual polyhedron.
// Erases v and re-triangulates the resulting hole using the recorded path.
template <class Polyhedron>
void erase_vertex_and_clip_ears(
  Polyhedron & dual,
  typename Polyhedron::Vertex_handle & v,
  const size_t start_vertex_id,
  const std::vector<int> & path);


// Collect all neighbor vertex handles of v in one-ring order.
template <class Vertex_handle>
std::vector<Vertex_handle>
collect_neighbors(const Vertex_handle & v);


// Compute the cost (primal volume added) of erasing dual vertex v, and record
// the ear-clipping path needed to re-triangulate its hole.
template <class Polyhedron>
std::pair<typename Polyhedron::Traits::FT, Record>
measure_vertex_erasure(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  const int max_degree_for_flips,
  const CostFunction cost_function = CostFunction::PRIMAL_VOLUME);

