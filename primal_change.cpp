#include "primal_change.hpp"
#include "primal_area_change.hpp"
#include "primal_volume_subtended.hpp"
// So we can static_assert that we're not using it.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, bool> primal_change(
  Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0,
  const CostFunction cost_function)
{
  using Scalar = typename Polyhedron::Traits::FT;
  using Point = typename Polyhedron::Point_3;

  ///////////////////////////////////////////////////////////////////
  std::vector<Point> all_primal_vertices;
  std::vector<bool> all_already;
  int num_boundary_halfedges = 0;
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
    const auto h0 = h;
    do
    {
      h->id() = id++;
      num_boundary_halfedges++;
      h = h->next();
    }while(h != h0);
    all_primal_vertices.resize(id);
    all_already.resize(id, false);
  }

  const auto dual_triangle_to_point = [&all_primal_vertices,&all_already](
      const int id,
      const Point & A, const Point & B, const Point & C)
  {
    if(all_already[id])
    {
      return all_primal_vertices[id];
    }
    // N = (A-B)×(C-B)
    const auto n = CGAL::cross_product(A - B, C - B);
    const auto bc = CGAL::centroid(A, B, C) - CGAL::ORIGIN;
    const auto beta = CGAL::scalar_product(n, bc);
    // n / beta
    const auto p = CGAL::ORIGIN + (n / beta);
    all_primal_vertices[id] = p;
    all_already[id] = true;
    return p;
  };

  const auto polygon_vector_area = [&all_primal_vertices](const std::vector<int> & vertices)
  {
    using Vector = typename CGAL::Kernel_traits<Point>::Kernel::Vector_3;

    Vector sum = CGAL::NULL_VECTOR;

    for (int i = 0; i < vertices.size(); i++)
    {
      const auto & v0 = all_primal_vertices[vertices[i]];
      const auto & v1 = all_primal_vertices[vertices[(i + 1) % vertices.size()]];

      sum = sum + CGAL::cross_product(v0 - CGAL::ORIGIN, v1 - CGAL::ORIGIN);
    }

    return sum;
  };
  const auto centroid = [&all_primal_vertices](const std::vector<int> & vertices)
  { 
    Point bc(0,0,0);
    for(const auto & v : vertices)
      bc = bc + (all_primal_vertices[v] - CGAL::ORIGIN);
    return ((bc-CGAL::ORIGIN) / Scalar(vertices.size()));
  };

  Scalar contribution_bottom = 0;

  for(typename Polyhedron::Vertex_const_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
  {
    //printf("v->id() = %d\n", v->id());
    const auto hv = v->halfedge();
    auto h = hv;
    std::vector<int> primal_vertices;primal_vertices.reserve(16);
    // consider each facet incident on the vertex
    do{
      //printf("  h: %d,%d | %s\n",
      // h->vertex()->id(),
      // h->opposite()->vertex()->id(),
      // h->is_border() ? "border" : "internal");
      if(h->is_border()) 
      {
        const auto q1 = 
          dual_triangle_to_point(
              h->id(),
              h->vertex()->point(),
              p0,
              h->opposite()->vertex()->point());
        //printf("%d: %g,%g,%g\n",h->id(), q1.x(), q1.y(), q1.z());
        primal_vertices.push_back(h->id());
        const auto q2 = 
          dual_triangle_to_point(
              h->next()->id(),
              h->vertex()->point(),
              h->next()->vertex()->point(),
              p0);
        //printf("%d: %g,%g,%g\n",h->next()->id(), q2.x(), q2.y(), q2.z());
        primal_vertices.push_back(h->next()->id());
      }else
      {
        const auto q = 
          dual_triangle_to_point(
              h->facet()->id(),
              h->vertex()->point(),
              h->next()->vertex()->point(),
              h->next()->next()->vertex()->point());
        //printf("%d: %g,%g,%g\n",h->facet()->id(), q.x(), q.y(), q.z());
        primal_vertices.push_back(h->facet()->id());
      }
      h = h->next()->opposite();
    } while(h != hv);
    // area of the primal polygon
    const auto n_v = polygon_vector_area(primal_vertices);
    switch(cost_function)
    {
      case CostFunction::PRIMAL_VOLUME:
      {
        const auto bc_v = centroid(primal_vertices);
        // volume += n⋅bc
        contribution_bottom += CGAL::scalar_product(n_v, bc_v);
        break;
      }
      case CostFunction::PRIMAL_AREA_CHANGE:
      {
        static_assert(!std::is_same_v<Scalar, CGAL::Epeck::FT>, "Scalar must not be Epeck");
        contribution_bottom += CGAL::sqrt(n_v.squared_length()) * 0.5;
        break;
      }
    }
  }

  Scalar contribution_top = 0;
  //printf("----------------------------------\n");
  {
    std::vector<int> primal_vertices;primal_vertices.reserve(num_boundary_halfedges);
    typename Polyhedron::Halfedge_const_iterator h = poly.halfedges_begin();
    for(; h != poly.halfedges_end(); ++h)
    {
      if(h->is_border()) break;
    }
    //printf("h: %d,%d\n", h->vertex()->id(), h->opposite()->vertex()->id());
    const auto h0 = h;
    do
    {
      //printf("h: %d,%d | %s\n",
      //  h->vertex()->id(),
      //  h->opposite()->vertex()->id(),
      //  h->is_border() ? "border" : "internal");
      const auto p = dual_triangle_to_point(
          h->id(),
          h->vertex()->point(),
          h->opposite()->vertex()->point(),
          p0);
      //printf("  %g, %g, %g\n", p.x(), p.y(), p.z());
      //printf("%d: %g,%g,%g\n",h->id(), p.x(), p.y(), p.z());
      primal_vertices.push_back(h->id());
      h = h->next();
    }while(h != h0);
    const auto n = polygon_vector_area(primal_vertices);
    switch(cost_function)
    {
      case CostFunction::PRIMAL_VOLUME:
      {
        const auto bc = centroid(primal_vertices);
        contribution_top += CGAL::scalar_product(n, bc);
        break;
      }
      case CostFunction::PRIMAL_AREA_CHANGE:
      {
        contribution_top += CGAL::sqrt(n.squared_length()) * 0.5;
        break;
      }
    }
  }
  Scalar cost;
  switch(cost_function)
  {
    case CostFunction::PRIMAL_VOLUME:
    {
      cost = (contribution_top - contribution_bottom)/Scalar(6);
      break;
    }
    case CostFunction::PRIMAL_AREA_CHANGE:
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
    const auto ori = CGAL::orientation(A, B, C, Point(0,0,0));
    if(ori != CGAL::ON_NEGATIVE_SIDE)
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
  
  // Could let the caller do this.
  cost = contains_origin ? std::numeric_limits<Scalar>::infinity() : cost ;
  return {cost, contains_origin};
}


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
template std::tuple<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, bool> primal_change<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>> &, CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Point_3 const&, CostFunction);
