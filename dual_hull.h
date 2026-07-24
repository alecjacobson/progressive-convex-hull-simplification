#pragma once
#include "geometry.h"
#include "polyhedron_utils.h"
#include "chebyshev_center.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>
#include <Eigen/Core>
#include <vector>
#include <tuple>

// Given a triangulated CGAL polyhedron, compute the Chebyshev center of the
// halfspaces corresponding to faces with squared area > primal_squared_area_tol.
// Converts to double to call igl::chebyshev_center, returns result in Kernel.
template <class Polyhedron>
typename Polyhedron::Traits::Point_3 polyhedron_chebyshev_center(
  const Polyhedron & poly,
  const double primal_squared_area_tol)
{
  using Kernel = typename Polyhedron::Traits;
  using Point = typename Kernel::Point_3;

  int nf = 0;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    if(geom::to_double(geom::squared_area(
        h->vertex()->point(),
        h->next()->vertex()->point(),
        h->next()->next()->vertex()->point())) > primal_squared_area_tol) { nf++; }
  }
  Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor> P(nf, 4);
  int i = 0;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const auto & A = h->vertex()->point();
    const auto & B = h->next()->vertex()->point();
    const auto & C = h->next()->next()->vertex()->point();
    if(geom::to_double(geom::squared_area(A,B,C)) <= primal_squared_area_tol) { continue; }
    const Eigen::Vector3d a(geom::to_double(A.x()),geom::to_double(A.y()),geom::to_double(A.z()));
    const Eigen::Vector3d b(geom::to_double(B.x()),geom::to_double(B.y()),geom::to_double(B.z()));
    const Eigen::Vector3d c(geom::to_double(C.x()),geom::to_double(C.y()),geom::to_double(C.z()));
    const Eigen::Vector3d n = (b-a).cross(c-a);
    const Eigen::Vector3d bc = (a+b+c)/3.0;
    P.row(i++) << n(0), n(1), n(2), -n.dot(bc);
  }
  Eigen::Matrix<double,3,1> x0d;
  igl::chebyshev_center(P, x0d);
  return Point(x0d(0), x0d(1), x0d(2));
}

// Given a triangulated CGAL polyhedron and interior point x0, return one dual
// point per face (skipping faces with squared area <= primal_squared_area_tol).
// Duality map: given halfspace [n, b] (b = -n·bc), d = -n / (n·x0 + b).
// All arithmetic stays in Kernel.
template <class Polyhedron>
std::vector<typename Polyhedron::Traits::Point_3> dual_points_list(
  const Polyhedron & poly,
  const double primal_squared_area_tol,
  const typename Polyhedron::Traits::Point_3 & x0)
{
  using Kernel = typename Polyhedron::Traits;
  using FT = typename Kernel::FT;
  using Point = typename Kernel::Point_3;
  using Vector = typename Kernel::Vector_3;

  std::vector<Point> dpts;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const auto & A = h->vertex()->point();
    const auto & B = h->next()->vertex()->point();
    const auto & C = h->next()->next()->vertex()->point();
    if(geom::to_double(geom::squared_area(A,B,C)) <= primal_squared_area_tol) { continue; }
    const Vector n = geom::cross_product(B - A, C - A);
    const Vector bc = ((A - geom::ORIGIN) + (B - geom::ORIGIN) + (C - geom::ORIGIN)) / FT(3);
    const FT b = -(n * bc);
    const FT denom = n * (x0 - geom::ORIGIN) + b;
    dpts.push_back(geom::ORIGIN + (-(n / denom)));
  }
  return dpts;
}

template <class Polyhedron, typename x0_type>
Eigen::Matrix<typename Polyhedron::Traits::FT, Eigen::Dynamic, 3, Eigen::RowMajor>
dual_to_primal_points(
  const Polyhedron & dual,
  const x0_type & x0_exact)
{
  using Scalar = typename Polyhedron::Traits::FT;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> pV(dual.size_of_facets(), 3);
  {
    for(auto f = dual.facets_begin(); f != dual.facets_end(); ++f)
    {
      auto h = f->halfedge();
      const auto & A = h->vertex()->point();
      const auto & B = h->next()->vertex()->point();
      const auto & C = h->next()->next()->vertex()->point();
      const auto n = geom::cross_product(B - A, C - A);
      const auto bc = ((A - geom::ORIGIN) + (B - geom::ORIGIN) + (C - geom::ORIGIN)) / Scalar(3);
      const Scalar beta = n * bc;
#if defined(PCHS_BACKEND_NATIVE)
      // Single kernel: x0 is already in the dual's kernel, no conversion.
      const typename Polyhedron::Traits::Point_3 & x0 = x0_exact;
#else
      using EK = typename CGAL::Kernel_traits<x0_type>::Kernel;
      typename Polyhedron::Traits::Point_3 x0;
      {
        CGAL::Cartesian_converter<EK, typename Polyhedron::Traits> to_dual_kernel;
        x0 = to_dual_kernel(x0_exact);
      }
#endif
      const auto p = geom::ORIGIN + (x0 - geom::ORIGIN) + (n / beta);
      pV.row(f->id()) << p.x(), p.y(), p.z();
    }
  }
  return pV;
}

// Convert a simplified dual polyhedron back to a primal polygon mesh.
// x0_exact is the Chebyshev center (may be in a different kernel than dual).
// Requires face ids to be set on dual (0..size_of_facets-1).
// Returns (pV, pPI, pPC): primal vertices, polygon index list, polygon count list.
template <class Polyhedron, typename x0_type>
std::tuple<
    Eigen::Matrix<typename Polyhedron::Traits::FT, Eigen::Dynamic, 3, Eigen::RowMajor>,
    Eigen::VectorXi,
    Eigen::VectorXi>
dual_to_primal_mesh(
  const Polyhedron & dual,
  const x0_type & x0_exact)
{
  const auto pV = dual_to_primal_points(dual, x0_exact);
  std::vector<int> pPI;
  std::vector<int> pPC = {0};
  {
    for(auto v = dual.vertices_begin(); v != dual.vertices_end(); ++v)
    {
      int np = 0;
      if(v->halfedge() == nullptr)
      {
#ifndef NDEBUG
        printf("Warning: vertex %d has no halfedge\n", v->id());
#endif
      }else
      {
        auto h_start = v->halfedge()->opposite()->prev();
        auto h = h_start;
        do
        {
          pPI.push_back(h->facet()->id());
          h = h->opposite()->prev();
          np++;
        } while(h != h_start);
      }
      pPC.push_back(pPC.back() + np);
    }
  }
  return {
    pV,
    Eigen::VectorXi(Eigen::VectorXi::Map(pPI.data(), pPI.size())),
    Eigen::VectorXi(Eigen::VectorXi::Map(pPC.data(), pPC.size()))};
}
