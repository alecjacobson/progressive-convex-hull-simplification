#include "primal_volume_subtended.hpp"
#include <igl/matlab_format.h>
#include <Eigen/Geometry>
#include <iostream>
#include <array>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template 
  <typename DerivedV,
   typename DerivedPI,
   typename DerivedPC>
std::tuple<typename DerivedV::Scalar, typename DerivedV::Scalar,bool> primal_volume_subtended(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedPI> & PI,
  const Eigen::MatrixBase<DerivedPC> & PC,
  const int i0)
{
  // Decomposition into tets using i0 as the origin seems numerically the best
  // bet.
  using Scalar = typename DerivedV::Scalar;
  Scalar primal_volume = 0;
  Scalar dual_volume = 0;
  // each polygonal face
  bool any_contains_origin = false;
  for(int p = 0;p<PC.size()-1;p++)
  {
    //printf("polygon %d\n",p);
    const int i = PC(p) + 0;
    // fan of triangles tessellating the polygonal face
    for(int j = PC(p) + 1;j<PC(p+1)-1;j++)
    {
      //printf("  triangle %d\n",j-PC(p));
      const int k = j + 1;
      // Tet with vertices i,j,k and i0. 
      using MatrixS43r = Eigen::Matrix<Scalar,4,3,Eigen::RowMajor>;
      const MatrixS43r Vpj =
        (MatrixS43r() << 
         V.row(PI(i)),
         V.row(PI(j)),
         V.row(PI(k)),
         V.row(i0)).finished();
      auto [pv,dv,contains_origin] = primal_volume_subtended(Vpj);
      if(contains_origin)
      {
        any_contains_origin = true;
        pv = std::numeric_limits<Scalar>::infinity();
      }else
      {
        const Scalar dv_tol = 1e-10;
        const Scalar pv_tol = 1e-15;
        // For convex simplification use case, I claim that ≈0 dual volume implies
        // primal volume should also be negligible. 
        //
        // It may be important to filter these per tet, in case some other
        // faces contribute non-negligible dual/primal volume.
        if(std::abs(dv) < dv_tol && std::abs(pv) > pv_tol)
        {
          pv = 0.0;
        }
      }
      primal_volume += pv;
      dual_volume += dv;
    }
  }
  return {primal_volume, dual_volume,any_contains_origin};
}

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, typename Polyhedron::Traits::FT, bool> primal_volume_subtended(
  const Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0)
{
  using Scalar = typename Polyhedron::Traits::FT;
  using MatrixS43r = Eigen::Matrix<Scalar,4,3,Eigen::RowMajor>;
  assert(!poly.is_closed() && "Should be a disk");

  Scalar dual_volume = 0;
  Scalar primal_volume = 0;
  bool any_contains_origin = false;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    const auto & A = f->halfedge()->vertex()->point();
    const auto & B = f->halfedge()->next()->vertex()->point();
    const auto & C = f->halfedge()->next()->next()->vertex()->point();
    const MatrixS43r Vpj = 
      (MatrixS43r() << 
       A.x(),A.y(),A.z(),
       B.x(),B.y(),B.z(),
       C.x(),C.y(),C.z(),
       p0.x(),p0.y(),p0.z()).finished();
    auto [pv,dv,contains_origin] = primal_volume_subtended(Vpj);
    if(contains_origin)
    {
      any_contains_origin = true;
      pv = std::numeric_limits<Scalar>::infinity();
    }else
    {
      const Scalar dv_tol = 1e-10;
      const Scalar pv_tol = 1e-15;
      // For convex simplification use case, I claim that ≈0 dual volume implies
      // primal volume should also be negligible. 
      //
      // It may be important to filter these per tet, in case some other
      // faces contribute non-negligible dual/primal volume.
      if(CGAL::abs(dv) < dv_tol && CGAL::abs(pv) > pv_tol)
      {
        pv = 0.0;
      }
    }
    primal_volume += pv;
    dual_volume += dv;
  }
  return {primal_volume, dual_volume,any_contains_origin};
}


template 
  <typename DerivedV>
std::tuple<typename DerivedV::Scalar, typename DerivedV::Scalar,bool> primal_volume_subtended(
  const Eigen::MatrixBase<DerivedV> & V)
{
  using Scalar = typename DerivedV::Scalar;
  // Compiler will unroll?
  static constexpr std::array<std::array<int,3>,4> F = {{
    {{0,1,2}},
    {{3,1,0}},
    {{3,0,2}},
    {{3,2,1}}
  }};
  Eigen::Matrix<Scalar,4,1> beta;
  using MatrixS33r = Eigen::Matrix<Scalar,3,3,Eigen::RowMajor>;
  Scalar beta_max = std::numeric_limits<Scalar>::lowest();
  for(int i = 0;i<4;i++)
  {
    const MatrixS33r M =
      (MatrixS33r() << 
       V.row(F[i][0]),
       V.row(F[i][1]),
       V.row(F[i][2])).finished();
    beta(i) = M.determinant();
    beta_max = std::max(beta_max,beta(i));
  }

  const MatrixS33r A =
    (MatrixS33r() << 
      V.row(0) - V.row(3),
      V.row(1) - V.row(3),
      V.row(2) - V.row(3)).finished();
  const Scalar delta = A.determinant();
  const Scalar dual_volume = delta / Scalar(6);

  const bool contains_origin = -beta_max > -1e-15;
  //printf("primal_volume_subtended_contains_origin = %d\n",contains_origin);
  Scalar primal_volume = 
      (delta * delta * delta) / (Scalar(6) * beta(0) * beta(1) * beta(2) * beta(3));
  //std::cout<<igl::matlab_format(beta,"beta")<<std::endl;
  //std::cout<<igl::matlab_format(delta,"delta")<<std::endl;
  //std::cout<<igl::matlab_format(primal_volume,"primal_volume")<<std::endl;
  //std::cout<<igl::matlab_format(dual_volume,"dual_volume")<<std::endl;

  //if(primal_volume < -1e-15)
  //{
  //  // Primal volume is suspiciously negative. 
  //  printf("     primal_volume = %g\n",primal_volume);
  //  printf("     delta = %g\n",delta);
  //  printf("     beta = [%g, %g, %g, %g]\n",beta(0),beta(1),beta(2),beta(3));
  //  printf("     dual_volume = %g\n",dual_volume);
  //  assert(std::abs(dual_volume) < 1e-15);
  //  // over-write with 0.0;
  //  primal_volume = 0.0;
  //}

  return {primal_volume, dual_volume,contains_origin};
}

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, bool> primal_volume_subtended_explicit(
  const Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0)
{
  //printf("primal_volume_subtended_explicit\n");
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

  const auto polygon_vector_area = [](const std::vector<Point> & vertices)
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
    return sum;
  };
  const auto centroid = [](const std::vector<Point> & vertices)
  { 
    Point bc(0,0,0);
    for(const auto & v : vertices)
      bc = bc + (v - CGAL::ORIGIN);
    return ((bc-CGAL::ORIGIN) / Scalar(vertices.size()));
  };
 
  Scalar volume_bottom = 0;
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
    const auto n_v = polygon_vector_area(primal_vertices);
    const auto bc_v = centroid(primal_vertices);
    // volume += n⋅bc
    volume_bottom += CGAL::scalar_product(n_v, bc_v);
  }
  Scalar volume_top = 0;
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
    const auto n = polygon_vector_area(primal_vertices);
    const auto bc = centroid(primal_vertices);
    volume_top += CGAL::scalar_product(n, bc);
  }
  const Scalar volume = (volume_top-volume_bottom)/Scalar(6);
  //printf("volume_bottom = %g\n", volume_bottom);
  //printf("volume_top = %g\n", volume_top);
  //printf("volume = %g\n", volume);

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

  return {volume , contains_origin};
}

template std::tuple<Eigen::Matrix<double, -1, 3, 1, -1, 3>::Scalar, Eigen::Matrix<double, -1, 3, 1, -1, 3>::Scalar, bool> primal_volume_subtended<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const&, int);
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
template std::tuple<
  CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, 
  CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, 
  bool> primal_volume_subtended<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>> const&, CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Point_3 const&);
template std::tuple<
  CGAL::Polyhedron_3<CGAL::Epeck, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, 
  CGAL::Polyhedron_3<CGAL::Epeck, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, 
  bool> primal_volume_subtended<CGAL::Polyhedron_3<CGAL::Epeck, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epeck, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>> const&, CGAL::Polyhedron_3<CGAL::Epeck, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Point_3 const&);
template std::tuple<
  CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT,
  bool> primal_volume_subtended_explicit<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>> const&, CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Point_3 const&);
