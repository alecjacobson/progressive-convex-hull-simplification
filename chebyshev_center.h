#ifndef IGL_CHEBYSHEV_CENTER_H
#define IGL_CHEBYSHEV_CENTER_H
#include <Eigen/Core>

namespace igl
{
  template <
   typename DerivedP,
   typename Derivedx0>
  bool chebyshev_center(
    const Eigen::MatrixBase<DerivedP>& P,
    Eigen::PlainObjectBase<Derivedx0>& x0);
}
#endif
