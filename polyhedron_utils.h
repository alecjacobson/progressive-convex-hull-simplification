#pragma once
#include <Eigen/Core>
#include <vector>

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
