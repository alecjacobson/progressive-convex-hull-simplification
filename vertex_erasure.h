#pragma once
#include "convex_triangulation.h"
#include "primal_volume_subtended.hpp"

#include <cassert>
#include <stdexcept>
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
clip_ears(Polyhedron & one_ring_copy, typename Polyhedron::Halfedge_handle h0)
{
  const auto walk = [](const auto h)->auto
  {
    assert(h->opposite()->is_border());
    return h->opposite()->prev()->opposite();
  };

#ifndef NDEBUG
  {
    int count = 0;
    auto h = h0;
    while(true)
    {
      h = walk(h);
      count++;
      assert(count < 1000);
      if(h == h0) break;
    }
    assert(count == (int)one_ring_copy.size_of_vertices());
  }
#endif

  int walks_since_last_ear = 0;
  auto h = h0;
  auto start_vertex_id = h->vertex()->id();
  std::vector<int> path;
  path.reserve(one_ring_copy.size_of_facets()-1);

  while(one_ring_copy.size_of_vertices() > 3)
  {
    auto w = walk(h);
    auto n = h->next();
    if(w == n)
    {
      // Ear found
      path.push_back(walks_since_last_ear);
      walks_since_last_ear = 0;
      auto q = h->next()->next()->opposite();
      one_ring_copy.erase_facet(h);
      assert(q->opposite()->is_border());
      h = q;
    }
    else
    {
      h = w;
      walks_since_last_ear++;
    }
  }
  return {start_vertex_id, path};
}

// Apply a previously recorded ear-clipping path to the actual dual polyhedron.
// Erases v and re-triangulates the resulting hole using the recorded path.
template <class Polyhedron>
void erase_vertex_and_clip_ears(
  Polyhedron & dual,
  typename Polyhedron::Vertex_handle & v,
  const size_t start_vertex_id,
  const std::vector<int> & path)
{
  auto h0 = dual.erase_center_vertex(v->halfedge());
  while(h0->vertex()->id() != start_vertex_id)
    h0 = h0->next();

  const auto walk = [](const auto h)->auto { return h->next(); };

  auto h = h0;
  for(int i = 0; i < (int)path.size(); i++)
  {
    const int steps = path[i];
    for(int j = 0; j < steps; j++)
      h = walk(h);
    h = dual.split_facet(h->prev(), h->next());
  }
}

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

// Compute the cost (primal volume added) of erasing dual vertex v, and record
// the ear-clipping path needed to re-triangulate its hole.
template <class Polyhedron>
std::pair<typename Polyhedron::Traits::FT, Record>
measure_vertex_erasure(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  const int max_degree_for_flips)
{
  using Scalar = typename Polyhedron::Traits::FT;
#ifndef NDEBUG
  {
    for(auto u = dual.vertices_begin(); u != dual.vertices_end(); ++u)
    {
      if(u->id() == v->id())
      {
        assert(u == v && "vertex ids should be unique");
        break;
      }
    }
  }
#endif

  auto p = v->point();
  Polyhedron one_ring_copy;
  typename Polyhedron::Halfedge_handle h0;

  if((int)v->degree() <= max_degree_for_flips)
    one_ring_triangulation_convex_via_flips(dual,v,one_ring_copy,h0);
  else
    one_ring_triangulation_convex_via_convex_hull(dual,v,one_ring_copy,h0);

#ifndef NDEBUG
  for(auto f = one_ring_copy.facets_begin(); f != one_ring_copy.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const auto ori = CGAL::orientation(
      h->vertex()->point(),
      h->next()->vertex()->point(),
      h->next()->next()->vertex()->point(),
      p);
    if(ori == CGAL::NEGATIVE)
      throw std::runtime_error("Removed point lies on negative side of a new face.");
  }
#endif

  auto [primal_volume, dual_volume, contains_origin] = primal_volume_subtended(one_ring_copy,p);
  const Scalar cost = primal_volume;

  Record record;
  std::tie(record.start_vertex_id, record.path) = clip_ears(one_ring_copy, h0);

  return {cost, record};
}
