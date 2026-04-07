#include "chebyshev_center.h"
#include <igl/PlainMatrix.h>
#include <igl/PlainVector.h>
#include <igl/matlab_format.h>
#define IGL_LINPROG_VERBOSE
#include <igl/linprog.h>
#include <sdlp/sdlp.hpp>
#include <iostream>

template <
 typename DerivedP,
 typename Derivedx0>
bool igl::chebyshev_center(
  const Eigen::MatrixBase<DerivedP>& P,
  Eigen::PlainObjectBase<Derivedx0>& x0)
{
  constexpr int N = DerivedP::ColsAtCompileTime;
  constexpr int M = DerivedP::RowsAtCompileTime;
  static_assert(N != Eigen::Dynamic, "DerivedP must have a fixed number of columns");

  // norm_vector = sqrt(sum(P(:,1:end-1).^2,2));
  using Scalar = typename DerivedP::Scalar;
  // static assert that Scalar is double
  static_assert(std::is_same<Scalar,double>::value,"Scalar must be double");
  igl::PlainVector<DerivedP,DerivedP::RowsAtCompileTime> norm_vector = P.leftCols(P.cols()-1).rowwise().norm();
  // Annoyingly igl::linprog is not templated well so use colmajor and double for all of
  // these.
  //
  // sdlp::linprog is also expecting colmajor...
  Eigen::Matrix<Scalar,N,1> c(N);
  c.setZero();
  c(c.size()-1) = -1;

  Eigen::Matrix<Scalar,M,N> A(P.rows(),N);
  A.leftCols(N-1) = P.leftCols(N-1);
  A.col(N-1) = norm_vector;
  Eigen::Matrix<Scalar,M,1> b = -P.col(N-1);
#define IGL_CHEBYSHEV_CENTER_USE_SDLP
#ifdef IGL_CHEBYSHEV_CENTER_USE_SDLP
  Eigen::Matrix<Scalar , N, 1> res;
  //std::cout<<igl::matlab_format(Eigen::MatrixXd(c),"c")<<std::endl;
  //std::cout<<igl::matlab_format(Eigen::MatrixXd(A),"A")<<std::endl;
  //std::cout<<igl::matlab_format(Eigen::MatrixXd(b),"b")<<std::endl;

  sdlp::linprog<N>(c,A,b,res);
#else
  assert(false && "This enforces x≥0 which we don't want...");
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> res;
  bool flag = igl::linprog( c,A,b,A.rows(),res);
  if(!flag) { return false; }
#endif
  Scalar constraint_tol = 1e-7;
  for(int i = 0;i<A.rows();i++)
  {
    if((A.row(i).dot(res)-b(i))>constraint_tol)
    {
      return false;
    }
  }

  x0 = res.template head<N-1>(N-1);
  return true;
}

template bool igl::chebyshev_center<Eigen::Matrix<double, -1, 4, 1, -1, 4>, Eigen::Matrix<double, 1, 3, 1, 1, 3>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 4, 1, -1, 4>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3>>&);
template bool igl::chebyshev_center<Eigen::Matrix<double, -1, 4, 1, -1, 4>, Eigen::Matrix<double, 3, 1, 0, 3, 1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, 4, 1, -1, 4>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1>>&);
