#pragma once
#include <Eigen/Core>
#include <igl/matlab_format.h>
#include <vector>
#include <iostream>

// Convert an Eigen matrix of 3D points to a vector of CGAL points.
template <typename Kernel, typename DerivedV>
std::vector<typename Kernel::Point_3> point_list(
    const Eigen::MatrixBase<DerivedV> & V)
{
  using Point = typename Kernel::Point_3;
  std::vector<Point> points(V.rows());
  for(int i = 0;i<V.rows();i++)
  {
    points[i] = Point(
      static_cast<typename Kernel::FT>(V(i,0)),
      static_cast<typename Kernel::FT>(V(i,1)),
      static_cast<typename Kernel::FT>(V(i,2))
      );
  }
  return points;
}

// Convert a vector of CGAL points to an Eigen matrix.
template <typename Kernel>
Eigen::Matrix<typename Kernel::FT,Eigen::Dynamic,3,Eigen::RowMajor> point_matrix(
    const std::vector<typename Kernel::Point_3> & points)
{
  Eigen::Matrix<typename Kernel::FT,Eigen::Dynamic,3,Eigen::RowMajor> V(points.size(),3);
  for(int i = 0;i<(int)points.size();i++)
  {
    V(i,0) = points[i].x();
    V(i,1) = points[i].y();
    V(i,2) = points[i].z();
  }
  return V;
}

template <class Polyhedron>
void print_VF(const Polyhedron & poly, const std::string prefix = "")
{
  int max_id = -1;
  //printf("print_VF\n");
  for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
  {
    //printf("v->id() = %zu\n",v->id());
    if(int(v->id()) > max_id) { max_id = int(v->id()); }
  }
  //printf("-----------------\n");
  //printf("max_id = %d\n",max_id);
  Eigen::Matrix<typename Polyhedron::Traits::FT,Eigen::Dynamic,3,Eigen::RowMajor> V(max_id+1,3);
  V.setConstant(std::numeric_limits<typename Polyhedron::Traits::FT>::quiet_NaN());
  for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    V.row(v->id()) << v->point().x(), v->point().y(), v->point().z();
  assert(poly.is_pure_triangle());
  Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> F(poly.size_of_facets(),3);
  {
    int k = 0;
    for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
    {
      auto h = f->halfedge();
      F.row(k++) <<
        h->vertex()->id(),
        h->next()->vertex()->id(),
        h->next()->next()->vertex()->id();
    }
    assert(k == (int)poly.size_of_facets());
  }
  std::cout<<igl::matlab_format(V,prefix+"V")<<std::endl;
  std::cout<<igl::matlab_format_index(F,prefix+"F")<<std::endl;
}
