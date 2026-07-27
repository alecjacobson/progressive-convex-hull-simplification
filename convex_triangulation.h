#pragma once
#include "geometry.h"
#include <cassert>
#include <functional>
#include <map>
#include <stdexcept>
#include <vector>

// Returns whether the shared edge between the two faces incident to halfedge e
// is convex (negative side), non-convex (positive side), or flat (boundary).
// Uses only exact orientation predicates — no plane constructions.
template <class Polyhedron>
typename Polyhedron::Traits::Oriented_side
edge_side(typename Polyhedron::Halfedge_const_handle e)
{
  using Point = typename Polyhedron::Point_3;
  using Halfedge_handle = typename Polyhedron::Halfedge_const_handle;

  Halfedge_handle h1 = e;
  Halfedge_handle h2 = e->opposite();

  const Point& a1 = h1->vertex()->point();
  const Point& b1 = h1->next()->vertex()->point();
  const Point& c1 = h1->next()->next()->vertex()->point();
  const Point& p2 = h2->next()->vertex()->point();

  const auto ori = geom::orientation(a1,b1,c1,p2);
  switch(ori)
  {
    case geom::POSITIVE: return geom::ON_POSITIVE_SIDE;
    case geom::NEGATIVE: return geom::ON_NEGATIVE_SIDE;
    default:             return geom::ON_ORIENTED_BOUNDARY;
  }
}

// Debug assertion: throws if any interior edge of poly is non-convex.
template <typename Polyhedron>
void confirm_all_edges_are_convex(const Polyhedron& poly)
{
  for(auto e = poly.halfedges_begin(); e != poly.halfedges_end(); ++e)
  {
    if(e > e->opposite()) continue;
    if(e->is_border_edge()) continue;
    if(edge_side<Polyhedron>(e) == geom::ON_POSITIVE_SIDE)
      throw std::runtime_error("Found a non-convex edge.");
  }
}

// Flip non-convex interior edges of poly until all are convex.
// Uses a stack to avoid re-scanning all edges after each flip.
// Returns number of flips performed. Throws if max_flips is exceeded.
template <class Polyhedron>
int flip_until_all_interior_edges_are_convex(
  Polyhedron& poly,
  const int max_flips = 1000)
{
  int num_flips = 0;
  std::vector<typename Polyhedron::Halfedge_handle> stack;
  stack.reserve(poly.size_of_halfedges()/2);

  // Push the canonical (smaller-index) half-edge of an interior edge, marked
  // "unresolved" (id = 0). Not a std::function -> no heap alloc / indirect call.
  auto mark_and_push_edge = [&stack](auto e)
  {
    if(e->opposite() < e) e = e->opposite();
    if(e->is_border_edge()) return;
    e->id() = false;
    stack.push_back(e);
  };

  for(auto e = poly.halfedges_begin(); e != poly.halfedges_end(); ++e)
  {
    if(e > e->opposite()) continue;
    mark_and_push_edge(e);
  }

  while(!stack.empty())
  {
    auto e = stack.back();
    stack.pop_back();
    if(e->id()) { continue; } // marked convex in the meantime
    auto side = edge_side<Polyhedron>(e);
    if(side == geom::ON_POSITIVE_SIDE)
    {
      mesh::flip_edge(e, poly);
      num_flips++;
      if(num_flips > max_flips)
        throw std::runtime_error("Too many flips, something is wrong.");
      mark_and_push_edge(e->next());
      mark_and_push_edge(e->next()->next());
      mark_and_push_edge(e->opposite()->next());
      mark_and_push_edge(e->opposite()->next()->next());
    }
  }

#ifndef NDEBUG
  confirm_all_edges_are_convex(poly);
#endif

  return num_flips;
}

// The remaining one-ring construction helpers use CGAL incremental-builder /
// make_triangle / convex_hull_3 machinery and are CGAL-backend only. The native
// backend builds one-rings via nat::Mesh::build + retriangulate_star instead
// (see native_algo.h), so these are not compiled or instantiated there.
#if !defined(PCHS_BACKEND_NATIVE)
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/convex_hull_3.h>

// Build a copy of the one-ring around vertex v (center + neighbors as a fan
// of triangles). Copies vertex ids. Returns by reference so that h0 remains
// a valid handle into one_ring_copy.
template <class Polyhedron>
Polyhedron extract_copy_of_one_ring(
  const Polyhedron& poly,
  typename Polyhedron::Vertex_const_handle v)
{
  using Point = typename Polyhedron::Point_3;
  Polyhedron out;

  auto h0 = v->halfedge();
  if(h0 == nullptr) return out;

  std::vector<Point> points;
  std::vector<std::size_t> vertex_ids;

  points.push_back(v->point());
  vertex_ids.push_back(v->id());

  auto h = h0;
  do
  {
    points.push_back(h->opposite()->vertex()->point());
    vertex_ids.push_back(h->opposite()->vertex()->id());
    h = h->next()->opposite();
  } while(h != h0);

  struct Builder : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
  {
    const std::vector<Point>& points;
    Builder(const std::vector<Point>& pts) : points(pts) {}

    void operator()(typename Polyhedron::HalfedgeDS& hds)
    {
      CGAL::Polyhedron_incremental_builder_3<typename Polyhedron::HalfedgeDS> B(hds, true);
      const int n = static_cast<int>(points.size()) - 1;
      B.begin_surface(points.size(), n);
      for(const auto& p : points) B.add_vertex(p);
      for(int i = 1; i <= n; ++i)
      {
        const int j = (i < n ? i+1 : 1);
        B.begin_facet();
        B.add_vertex_to_facet(0);
        B.add_vertex_to_facet(j);
        B.add_vertex_to_facet(i);
        B.end_facet();
      }
      B.end_surface();
    }
  };

  Builder builder(points);
  out.delegate(builder);

  {
    auto vv = out.vertices_begin();
    for(std::size_t i = 0; i < vertex_ids.size(); ++i, ++vv)
      vv->id() = vertex_ids[i];
  }
  return out;
}

// Fan-triangulate the one-ring of v into one_ring_copy, and set h0 to a
// canonical border halfedge for use with clip_ears / erase_vertex_and_clip_ears.
template <class Polyhedron>
void one_ring_triangulation(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy,
  typename Polyhedron::Halfedge_handle & h0)
{
  using Point = typename Polyhedron::Point_3;
  const int nv = v->degree();
  one_ring_copy.reserve(nv, 2*nv + 2*(nv-3), nv-2);
  // CGAL crashing on ->opposite() here sometimes...
  // Hard to reproduce in debug mode: Ta Prohm.stl, Processing Pan.stl
  const auto h_start = v->halfedge()->opposite()->prev();
  auto h = h_start;
  int i = 0;
  std::vector<typename Polyhedron::Vertex_const_handle> prev_verts(2);
  typename Polyhedron::Halfedge_handle q,H;
  do
  {
    const auto vert = h->opposite()->vertex();
    if(i < 2)
    {
      prev_verts[i] = vert;
    }
    else if(i == 2)
    {
      q = one_ring_copy.make_triangle(
        prev_verts[0]->point(),
        prev_verts[1]->point(),
                 vert->point());
      q->vertex()->id() = prev_verts[0]->id();
      q->next()->vertex()->id() = prev_verts[1]->id();
      q->next()->next()->vertex()->id() = vert->id();
      H = q->next()->opposite();
      assert(H->vertex() == q->vertex());
      assert(H->is_border());
    }
    else if(i > 2)
    {
      //   0------3
      //  ↑| ↖\   ↑
      // H||  qG  s
      //  |↓   \↘ |
      //   1----→2
      auto G = q->opposite();
      assert(G->is_border());
      auto s = one_ring_copy.add_vertex_and_facet_to_border(H,G);
      assert(s->opposite()->is_border());
      s->vertex()->point() = vert->point();
      s->vertex()->id() = vert->id();
      q = s->next();
    }
    i++;
    h = h->opposite()->prev();
  } while(h != h_start);
  h0 = q;
}

// Triangulate the one-ring of v convexly using edge flips.
template <class Polyhedron>
void one_ring_triangulation_convex_via_flips(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy,
  typename Polyhedron::Halfedge_handle & h0)
{
  one_ring_triangulation(dual,v,one_ring_copy,h0);
  const int num_flips = flip_until_all_interior_edges_are_convex(one_ring_copy,1000);
  assert(num_flips < 1000);
#ifndef NDEBUG
  confirm_all_edges_are_convex(one_ring_copy);
#endif
}

// Triangulate the one-ring of v convexly via convex hull of the neighbors.
// Used as fallback for high-degree vertices.
template <class Polyhedron>
void one_ring_triangulation_convex_via_convex_hull(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy,
  typename Polyhedron::Halfedge_handle & h0)
{
  using Point = typename Polyhedron::Point_3;
  auto p = v->point();

  std::map<Point,int> point_to_id;
  std::vector<Point> points;
  points.reserve(v->degree());
  {
    auto h = v->halfedge();
    do
    {
      auto point = h->opposite()->vertex()->point();
      auto id = h->opposite()->vertex()->id();
      point_to_id[point] = id;
      points.push_back(point);
      h = h->next()->opposite();
    } while(h != v->halfedge());
  }

  CGAL::convex_hull_3(points.begin(),points.end(),one_ring_copy);
  for(auto vv = one_ring_copy.vertices_begin(); vv != one_ring_copy.vertices_end(); ++vv)
    vv->id() = point_to_id[vv->point()];

  // Remove faces on the wrong side of the removed vertex p
  while(true)
  {
    bool found_any = false;
    for(auto f = one_ring_copy.facets_begin(); f != one_ring_copy.facets_end(); )
    {
      auto hh = f->halfedge();
      f++;
      const auto & a = hh->vertex()->point();
      const auto & b = hh->next()->vertex()->point();
      const auto & c = hh->next()->next()->vertex()->point();
      if(geom::orientation(a,b,c,p) == geom::NEGATIVE)
      {
        found_any = true;
        one_ring_copy.erase_facet(hh);
      }
    }
    if(!found_any) break;
  }

  // Find canonical h0 matching the first neighbor of v
  for(auto hh = one_ring_copy.halfedges_begin(); hh != one_ring_copy.halfedges_end(); ++hh)
  {
    if(hh->opposite()->is_border() &&
       hh->opposite()->vertex()->id() == v->halfedge()->opposite()->vertex()->id())
    {
      h0 = hh;
      break;
    }
  }
  assert(one_ring_copy.size_of_vertices() == v->degree());
}

#endif  // !PCHS_BACKEND_NATIVE
