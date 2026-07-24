#pragma once
#include "CostFunction.hpp"
#include "geometry.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

/// Cost of removing the vertex that used to correspond to p0, given `poly`: a
/// copy of the dual one-ring after that vertex was removed and the hole was
/// re-triangulated convexly. Returns (cost, contains_origin).
///
/// Header-only (templated on the mesh) so it can be instantiated on either
/// backend's polyhedron. All geometry goes through geom::; all mesh access is
/// the CGAL-Polyhedron-compatible navigation API.
template <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, bool> primal_change(
  Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0,
  const CostFunction cost_function = CostFunction::PRIMAL_VOLUME)
{
  using Scalar = typename Polyhedron::Traits::FT;
  using Point = typename Polyhedron::Point_3;

  ///////////////////////////////////////////////////////////////////
  std::vector<Point> all_primal_vertices;
  std::vector<bool> all_already;
  int num_boundary_halfedges = 0;
  typename Polyhedron::Halfedge_handle h0;
  {
    int id = 0;
    // First label all of the facets
    for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
    {
      f->id() = id++;
    }
    // Then label the boundary half-edges
    typename Polyhedron::Halfedge_iterator h = poly.halfedges_begin();
    for(; h != poly.halfedges_end(); ++h)
    {
      if(h->is_border()) break;
    }
    const auto h0_opposite = h;
    do
    {
      h->id() = id++;
      num_boundary_halfedges++;
      h = h->next();
    }while(h != h0_opposite);
    h0 = h0_opposite->opposite();
    all_primal_vertices.resize(id);
    all_already.resize(id, false);
  }

  const auto dual_triangle_to_point_uncached = [](
      const Point & A, const Point & B, const Point & C)
  {
    // N = (A-B)×(C-B)
    const auto n = geom::cross_product(A - B, C - B);
    const auto bc = geom::centroid(A, B, C) - geom::ORIGIN;
    const auto beta = geom::scalar_product(n, bc);
    // n / beta
    const auto p = geom::ORIGIN + (n / beta);
    return p;
  };

  const auto dual_triangle_to_point = [&all_primal_vertices,&all_already,&dual_triangle_to_point_uncached](
      const int id,
      const Point & A, const Point & B, const Point & C)
  {
    if(all_already[id])
    {
      return all_primal_vertices[id];
    }
    const auto p = dual_triangle_to_point_uncached(A,B,C);
    all_primal_vertices[id] = p;
    all_already[id] = true;
    return p;
  };

  const auto polygon_vector_area = [&all_primal_vertices](const std::vector<int> & vertices)
  {
    using Vector = typename Polyhedron::Traits::Vector_3;

    Vector sum = geom::NULL_VECTOR;

    for (int i = 0; i < vertices.size(); i++)
    {
      const auto & v0 = all_primal_vertices[vertices[i]];
      const auto & v1 = all_primal_vertices[vertices[(i + 1) % vertices.size()]];

      sum = sum + geom::cross_product(v0 - geom::ORIGIN, v1 - geom::ORIGIN);
    }

    return sum;
  };
  const auto centroid = [&all_primal_vertices](const std::vector<int> & vertices)
  {
    Point bc(0,0,0);
    for(const auto & v : vertices)
      bc = bc + (all_primal_vertices[v] - geom::ORIGIN);
    return ((bc-geom::ORIGIN) / Scalar(vertices.size()));
  };

  Scalar contribution_bottom = 0;

  for(typename Polyhedron::Vertex_const_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
  {
    const auto hv = v->halfedge();
    auto h = hv;
    std::vector<int> primal_vertices;primal_vertices.reserve(16);
    // consider each facet incident on the vertex
    do{
      if(h->is_border())
      {
        const auto q1 =
          dual_triangle_to_point(
              h->id(),
              h->vertex()->point(),
              p0,
              h->opposite()->vertex()->point());
        primal_vertices.push_back(h->id());
        const auto q2 =
          dual_triangle_to_point(
              h->next()->id(),
              h->vertex()->point(),
              h->next()->vertex()->point(),
              p0);
        primal_vertices.push_back(h->next()->id());
      }else
      {
        const auto q =
          dual_triangle_to_point(
              h->facet()->id(),
              h->vertex()->point(),
              h->next()->vertex()->point(),
              h->next()->next()->vertex()->point());
        primal_vertices.push_back(h->facet()->id());
      }
      h = h->next()->opposite();
    } while(h != hv);
    // area of the primal polygon
    const auto n_v = polygon_vector_area(primal_vertices);
    switch(cost_function)
    {
      case CostFunction::PRIMAL_MEAN_WIDTH: { /* handled elsewhere */ break; }
      case CostFunction::PRIMAL_VOLUME:
      {
        const auto bc_v = centroid(primal_vertices);
        // volume += n⋅bc
        contribution_bottom += geom::scalar_product(n_v, bc_v);
        break;
      }
      case CostFunction::PRIMAL_AREA:
      {
        static_assert(std::is_floating_point_v<Scalar>,
                      "area cost needs a floating-point (sqrt-able) scalar");
        contribution_bottom += geom::sqrt(n_v.squared_length()) * 0.5;
        break;
      }
    }
  }

  Scalar contribution_top = 0;
  {
    std::vector<int> primal_vertices;primal_vertices.reserve(num_boundary_halfedges);
    typename Polyhedron::Halfedge_const_iterator h = poly.halfedges_begin();
    for(; h != poly.halfedges_end(); ++h)
    {
      if(h->is_border()) break;
    }
    const auto h1 = h;
    do
    {
      const auto p = dual_triangle_to_point(
          h->id(),
          h->vertex()->point(),
          h->opposite()->vertex()->point(),
          p0);
      primal_vertices.push_back(h->id());
      h = h->next();
    }while(h != h1);
    const auto n = polygon_vector_area(primal_vertices);
    switch(cost_function)
    {
      case CostFunction::PRIMAL_MEAN_WIDTH: { /* handled elsewhere */ break; }
      case CostFunction::PRIMAL_VOLUME:
      {
        const auto bc = centroid(primal_vertices);
        contribution_top += geom::scalar_product(n, bc);
        break;
      }
      case CostFunction::PRIMAL_AREA:
      {
        contribution_top += geom::sqrt(n.squared_length()) * 0.5;
        break;
      }
    }
  }
  Scalar cost;
  switch(cost_function)
  {
    case CostFunction::PRIMAL_MEAN_WIDTH: { /* handled elsewhere */ break; }
    case CostFunction::PRIMAL_VOLUME:
    {
      cost = (contribution_top - contribution_bottom)/Scalar(6);
      break;
    }
    case CostFunction::PRIMAL_AREA:
    {
      cost = contribution_bottom - contribution_top;
      break;
    }
  }

  // Check that origin is on negative side of all faces
  bool contains_origin = false;
  for(typename Polyhedron::Facet_const_iterator f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    const auto & A = f->halfedge()->vertex()->point();
    const auto & B = f->halfedge()->next()->vertex()->point();
    const auto & C = f->halfedge()->next()->next()->vertex()->point();
    const auto ori = geom::orientation(A, B, C, Point(0,0,0));
    if(ori != geom::ON_NEGATIVE_SIDE)
    {
      contains_origin = true;
      break;
    }
  }

  if(cost_function == CostFunction::PRIMAL_VOLUME)
  {
    // numerically contains origin (otherwise cost will be nuts)
    contains_origin = contains_origin || cost < -1e-15;
    cost = (cost <= Scalar(0.0)) ? Scalar(0.0) : cost;
  }

  if(cost_function == CostFunction::PRIMAL_MEAN_WIDTH)
  {
    cost = 0;

    const auto mean_width_contribution = [](
        const auto & p_ni, const auto & p_nj, const auto & p_pu, const auto & p_pv)
    {
      // edge length between u and v
      const auto edge_vec = p_pu - p_pv;
      const auto edge_len = geom::sqrt(edge_vec.squared_length());
      if(edge_len == 0.0)
      {
        return Scalar(0.0);
      }
      const auto edge_dir = edge_vec / edge_len;
      //// dihedral angle between faces i and j
      const auto dihedral_angle = std::atan2(
          geom::scalar_product(geom::cross_product(p_nj, p_ni), edge_dir),
          geom::scalar_product(p_ni, p_nj));
      return edge_len * dihedral_angle;
    };

    const auto mean_width_contribution_internal_edge_to_p0 =
      [&mean_width_contribution,&dual_triangle_to_point_uncached,&p0](
        const auto & h)
    {
      //  a
      //  | \
      //  ↓h \
      //  |   \
      //  b----p₀
      //  |   /
      //  ↓h_next
      //  | /
      //  c
      const auto h_next = h->opposite()->prev()->opposite();
      const auto a = h->opposite()->vertex()->point();
      const auto b = h->vertex()->point();
      const auto c = h_next->vertex()->point();
      const auto p_ni = (p0 - geom::ORIGIN);
      const auto p_nj = (b - geom::ORIGIN);
      const auto p_pu = dual_triangle_to_point_uncached( a,b,p0);
      const auto p_pv = dual_triangle_to_point_uncached( b,c,p0);
      return mean_width_contribution(p_ni, p_nj, p_pu, p_pv);
    };

    const auto mean_width_contribution_internal_edge = [&mean_width_contribution,&dual_triangle_to_point_uncached](
        const auto & e)
    {
      assert(!e->is_border() && !e->opposite()->is_border() && "edge should be internal");
      const auto d_pi = e->vertex()->point();
      const auto d_pj = e->opposite()->vertex()->point();
      // primal face normals (unnormalized if we used atan2 in
      // mean_width_contribution)
      const auto p_ni = (d_pi - geom::ORIGIN);
      const auto p_nj = (d_pj - geom::ORIGIN);
      // (ij) between dual vertices i and j. Lies between dual faces u and v
      // which are always triangles and correspond to primal vertices u and v.
      const auto p_pu =
          dual_triangle_to_point_uncached(
              e->vertex()->point(),
              e->next()->vertex()->point(),
              e->next()->next()->vertex()->point());
      const auto p_pv =
          dual_triangle_to_point_uncached(
              e->opposite()->vertex()->point(),
              e->opposite()->next()->vertex()->point(),
              e->opposite()->next()->next()->vertex()->point());
      return mean_width_contribution(p_ni, p_nj, p_pu, p_pv);
    };

    for (auto e = poly.edges_begin(); e != poly.edges_end(); ++e)
    {
      if(e->is_border() || e->opposite()->is_border())
      {
        // skip border edges. We'll handle those separately.
        continue;
      }
      // internal edge from dual vertex i to j. This corresponds to primal faces
      // i and j. Between which there is a dihedral angle θij
      const auto contrib = mean_width_contribution_internal_edge(e);
      cost += contrib;
    }

    // "rim" contributions for both. "spoke" contributions for old.
    assert(h0->opposite()->is_border() && "h0 should be opposite a border half-edge");
    auto h = h0;
    do
    {
      const auto d_pi = h->vertex()->point();
      const auto d_pj = h->opposite()->vertex()->point();
      const auto p_ni = (d_pi - geom::ORIGIN);
      const auto p_nj = (d_pj - geom::ORIGIN);

      const auto p_pu_new =
          dual_triangle_to_point_uncached(
              h->vertex()->point(),
              h->next()->vertex()->point(),
              h->next()->next()->vertex()->point());
      const auto p_pu_old =
          dual_triangle_to_point_uncached(
              h->vertex()->point(),
              p0,
              h->next()->next()->vertex()->point());

      // We need to use the half-edge on the original polyhedron to get the
      // opposite facet
      assert(h->opposite()->is_border());
      const auto old_spoke = mean_width_contribution_internal_edge_to_p0(h);
      const auto flap_change = mean_width_contribution(p_ni, p_nj, p_pu_new, p_pu_old);
      cost += flap_change - old_spoke;
      h = h->opposite()->prev()->opposite();
    }while(h != h0);
    assert(h == h0);
    cost /= (4.0 * std::acos(-1.0));

    if(cost < Scalar(0.0))
    {
      assert(cost >= Scalar(-1));
    }
    cost = (cost <= Scalar(0.0)) ? Scalar(0.0) : cost;
  }

  // Could let the caller do this.
  cost = contains_origin ? std::numeric_limits<Scalar>::infinity() : cost ;
  return {cost, contains_origin};
}
