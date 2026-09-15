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

// Convenience wrapper around igl::chebyshev_center: given a set of halfspaces
// [nx,ny,nz,b] (interior: n·x + b <= 0), returns the center of the largest
// inscribed ball. Throws std::runtime_error if halfspaces isn't N x 4, or if
// the underlying LP fails outright. Note: the LP does not constrain the ball
// radius to be non-negative, so some genuinely infeasible halfspace sets
// (most commonly a pair of disjoint, opposing planes) can silently return a
// center with negative radius instead of throwing — the caller is
// responsible for sanity-checking the result if the input isn't already
// known to bound a nonempty region (e.g. a polytope's own faces).
Eigen::Vector3d chebyshev_center(const Eigen::MatrixXd & halfspaces);

#endif
