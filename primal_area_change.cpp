#include "primal_area_change.hpp"
#include <cstdio>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// So we can static_assert that we're not using it.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, typename Polyhedron::Traits::FT, bool> primal_area_change(
  const Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0)
{
  using Scalar = typename Polyhedron::Traits::FT;
  using Point = typename Polyhedron::Point_3;

  const auto dual_triangle_to_point = [](const Point & A, const Point & B, const Point & C)
  {
    // N = (A-B)×(C-B)
    const auto n = CGAL::cross_product(A - B, C - B);
    const auto bc = CGAL::centroid(A, B, C) - CGAL::ORIGIN;
    const auto beta = CGAL::scalar_product(n, bc);
    // n / beta
    return CGAL::ORIGIN + (n / beta);
  };

  const auto polygon_area = [](const std::vector<Point> & vertices)
  {
    using Vector = typename CGAL::Kernel_traits<Point>::Kernel::Vector_3;

    Vector sum = CGAL::NULL_VECTOR;

    for (int i = 0; i < vertices.size(); i++)
    {
      const auto & v0 = vertices[i];
      const auto & v1 = vertices[(i + 1) % vertices.size()];

      sum = sum + CGAL::cross_product(v0 - CGAL::ORIGIN, v1 - CGAL::ORIGIN);
    }

    // Static assert that Scalar is not Epeck
    static_assert(!std::is_same_v<Scalar, CGAL::Epeck::FT>, "Scalar must not be Epeck");
    return CGAL::sqrt(sum.squared_length()) * 0.5;
  };

  Scalar area_before = 0;
  Scalar area_after = 0;
  // consider each vertex
  for(typename Polyhedron::Vertex_const_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
  {
    //printf("v->id() = %d\n", v->id());
    const auto hv = v->halfedge();
    auto h = hv;
    std::vector<Point> primal_vertices;
    // consider each facet incident on the vertex
    do{
      //printf("  h: %d,%d | %s\n",
      // h->vertex()->id(),
      // h->opposite()->vertex()->id(),
      // h->is_border() ? "border" : "internal");
      if(h->is_border()) 
      {
        primal_vertices.push_back(dual_triangle_to_point(
            h->vertex()->point(),
            p0,
            h->opposite()->vertex()->point()));
        primal_vertices.push_back(dual_triangle_to_point(
            h->vertex()->point(),
            h->next()->vertex()->point(),
            p0));
      }else
      {
        primal_vertices.push_back(dual_triangle_to_point(
            h->vertex()->point(),
            h->next()->vertex()->point(),
            h->next()->next()->vertex()->point()));
      }
      h = h->next()->opposite();
    } while(h != hv);
    // area of the primal polygon
    Scalar area_v = polygon_area(primal_vertices);
    //printf("  area_v = %g\n", area_v);
    area_after += area_v;
  }
  // consider each half-edge
  {
    std::vector<Point> primal_vertices;
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
          h->vertex()->point(),
          h->opposite()->vertex()->point(),
          p0);
      //printf("  %g, %g, %g\n", p.x(), p.y(), p.z());
      primal_vertices.push_back(p);
      h = h->next();
    }while(h != h0);
    area_before = polygon_area(primal_vertices);
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

  return {area_before, area_after, contains_origin};
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

template std::tuple<
  CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, 
  CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, 
  bool> 
  primal_area_change<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>> const&, CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Point_3 const&);
