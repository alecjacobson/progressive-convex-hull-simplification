// Tier 1 — global invariant oracles. These must hold for ANY backend, at every
// step, with no golden files. They are the safety net for the refactor and for
// the eventual native backend.
#include "test_framework.h"
#include "test_backend.h"

#include <igl/icosahedron.h>
#include <igl/read_triangle_mesh.h>

#include <string>

namespace
{
// Euler characteristic of a closed triangle mesh (V - E + F), with E = 3F/2.
int euler_characteristic_closed_triangle(const Eigen::MatrixXd & V,
                                         const Eigen::MatrixXi & F)
{
  CHECK(F.cols() == 3);
  CHECK((3 * F.rows()) % 2 == 0);  // closed => even number of directed edges
  const int E = (3 * F.rows()) / 2;
  return int(V.rows()) - E + int(F.rows());
}

// For the primal polygon hull (pV,pPI,pPC), verify every input point in P lies
// inside (within tol) the halfspace of every face. This is the conservative /
// exterior-simplification guarantee: the hull strictly contains the input.
void check_contains_points(const Eigen::MatrixXd & pV,
                           const Eigen::VectorXi & pPI,
                           const Eigen::VectorXi & pPC,
                           const Eigen::MatrixXd & P,
                           double scale)
{
  // Interior reference point: centroid of hull vertices (interior for a convex
  // polytope with >0 volume).
  const Eigen::RowVector3d x0 = pV.colwise().mean();
  const double tol = 1e-6 * scale;

  for(int f = 0; f + 1 < pPC.size(); ++f)
  {
    const int start = pPC(f), n = pPC(f + 1) - pPC(f);
    CHECK(n >= 3);

    // Newell's method: robust polygon normal + a point on the plane.
    Eigen::RowVector3d nrm = Eigen::RowVector3d::Zero();
    Eigen::RowVector3d c   = Eigen::RowVector3d::Zero();
    for(int i = 0; i < n; ++i)
    {
      const Eigen::RowVector3d a = pV.row(pPI(start + i));
      const Eigen::RowVector3d b = pV.row(pPI(start + (i + 1) % n));
      nrm.x() += (a.y() - b.y()) * (a.z() + b.z());
      nrm.y() += (a.z() - b.z()) * (a.x() + b.x());
      nrm.z() += (a.x() - b.x()) * (a.y() + b.y());
      c += a;
    }
    c /= n;
    const double nlen = nrm.norm();
    if(nlen < 1e-20 * scale) continue;  // degenerate face, skip
    nrm /= nlen;
    // Orient outward: interior point must be on the negative side.
    if(nrm.dot(x0 - c) > 0) nrm = -nrm;
    const double b_off = nrm.dot(c);

    for(int p = 0; p < P.rows(); ++p)
    {
      const double s = nrm.dot(P.row(p)) - b_off;
      CHECK_MSG(s <= tol,
        "input point " + std::to_string(p) + " outside face " +
        std::to_string(f) + " by " + std::to_string(s));
    }
  }
}

double bbox_diag(const Eigen::MatrixXd & V)
{
  return (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
}

// Drive chs to a target and assert all structural invariants on the dual and the
// containment guarantee on the primal.
void check_invariants_at(ConvexHullSimplification & chs,
                         const Eigen::MatrixXd & inputV, int target)
{
  chs.simplify_to(target);

  auto [dV, dF] = chs.get_dual_mesh();
  CHECK_MSG(euler_characteristic_closed_triangle(dV, dF) == 2,
            "dual is not a closed genus-0 triangle mesh");
  CHECK(dV.rows() == chs.num_dual_vertices());
  CHECK(dF.rows() == 2 * dV.rows() - 4);  // closed triangulation

  auto [pV, pPI, pPC] = chs.get_primal_mesh();
  CHECK(pPC.size() - 1 == chs.num_dual_vertices());  // one polygon per dual vtx
  check_contains_points(pV, pPI, pPC, inputV, bbox_diag(inputV));
}
}  // namespace

TEST(invariants_icosahedron_each_step)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  ConvexHullSimplification chs(V, F);

  const int n0 = chs.num_dual_vertices();
  int prev = n0;
  // Every successful step removes exactly one dual vertex.
  while(chs.step())
  {
    const int now = chs.num_dual_vertices();
    CHECK_MSG(now == prev - 1, "a step did not remove exactly one vertex");
    prev = now;

    auto [dV, dF] = chs.get_dual_mesh();
    CHECK(euler_characteristic_closed_triangle(dV, dF) == 2);
    if(now <= 4) break;  // a tetrahedron is the minimal closed hull
  }
}

TEST(invariants_icosahedron_containment)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  for(int target : {18, 12, 8, 6, 4})
  {
    ConvexHullSimplification chs(V, F);
    check_invariants_at(chs, V, target);
  }
}

TEST(invariants_costs_nonnegative)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  ConvexHullSimplification chs(V, F, 100, CostFunction::PRIMAL_VOLUME);
  chs.simplify_to(6);
  const Eigen::VectorXd costs = chs.popped_dual_vertex_costs();
  const Eigen::VectorXi ids   = chs.popped_dual_vertex_ids();
  double cumulative = 0.0;
  for(int i = 0; i < ids.size(); ++i)
  {
    const double c = costs(ids(i));
    CHECK_MSG(c >= -1e-12, "popped cost is negative");
    cumulative += c;  // volume added is non-decreasing
  }
  CHECK(cumulative >= -1e-12);
}

// Actaeon is large; only run when the mesh is present next to the binary's cwd.
TEST(invariants_actaeon_if_present)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::read_triangle_mesh(std::string(PCHS_SOURCE_DIR) + "/Actaeon.ply", V, F))
  {
    // Not a failure — just nothing to check.
    return;
  }
  ConvexHullSimplification chs(V, F);
  check_invariants_at(chs, V, 200);
}
