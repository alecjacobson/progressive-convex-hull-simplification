// Phase B, brick 1 — validate the native geometry kernel (nat::) against CGAL's
// on exact-integer inputs (where both are exact, so signs/values must match
// bit-for-bit up to floating rounding of the constructions).
#include "test_framework.h"

#include "native_geom.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <array>
#include <vector>

namespace
{
using CK = CGAL::Exact_predicates_inexact_constructions_kernel;
using CP = CK::Point_3;
using CV = CK::Vector_3;

nat::Point3 np(double x, double y, double z) { return {x, y, z}; }
CP          cp(double x, double y, double z) { return CP(x, y, z); }

// A deterministic spread of small-integer points (exact as doubles).
std::vector<std::array<double, 3>> sample_points()
{
  std::vector<std::array<double, 3>> pts;
  for(int x = -3; x <= 3; ++x)
    for(int y = -3; y <= 3; ++y)
      for(int z = -3; z <= 3; ++z)
        pts.push_back({double(x), double(y), double(z)});
  return pts;
}
}  // namespace

TEST(native_geom_orient3d_matches_cgal)
{
  const auto pts = sample_points();
  // Fixed quadruples chosen by a deterministic stride so we cover positive,
  // negative, and coplanar (ZERO) configurations without RNG.
  const int n = int(pts.size());
  int checked = 0, zeros = 0, pos = 0, neg = 0;
  for(int i = 0; i < n; i += 7)
    for(int j = 1; j < n; j += 13)
      for(int k = 2; k < n; k += 17)
        for(int l = 3; l < n; l += 101)
        {
          const auto & A = pts[i]; const auto & B = pts[j];
          const auto & C = pts[k]; const auto & D = pts[l];
          const int cgal = CGAL::orientation(cp(A[0], A[1], A[2]),
                                             cp(B[0], B[1], B[2]),
                                             cp(C[0], C[1], C[2]),
                                             cp(D[0], D[1], D[2]));
          const int nate = nat::orient3d(np(A[0], A[1], A[2]),
                                         np(B[0], B[1], B[2]),
                                         np(C[0], C[1], C[2]),
                                         np(D[0], D[1], D[2]));
          CHECK(cgal == nate);
          checked++;
          if(nate == 0) zeros++; else if(nate > 0) pos++; else neg++;
        }
  CHECK(checked > 1000);
  CHECK(zeros > 0);  // exercised coplanar case
  CHECK(pos > 0 && neg > 0);
}

TEST(native_geom_constructions_match_cgal)
{
  const auto pts = sample_points();
  const int n = int(pts.size());
  const double tol = 1e-12;
  for(int i = 0; i < n; i += 11)
  {
    const auto & A = pts[i];
    const auto & B = pts[(i + 37) % n];
    const auto & C = pts[(i + 71) % n];

    const CP a = cp(A[0], A[1], A[2]), b = cp(B[0], B[1], B[2]), c = cp(C[0], C[1], C[2]);
    const nat::Point3 na = np(A[0], A[1], A[2]), nb = np(B[0], B[1], B[2]), nc = np(C[0], C[1], C[2]);

    // cross product of edge vectors
    const CV cx = CGAL::cross_product(b - a, c - a);
    const nat::Vector3 nx = nat::cross_product(nb - na, nc - na);
    CHECK_NEAR(CGAL::to_double(cx.x()), nx.x(), tol);
    CHECK_NEAR(CGAL::to_double(cx.y()), nx.y(), tol);
    CHECK_NEAR(CGAL::to_double(cx.z()), nx.z(), tol);

    // dot product (spelled operator* on vectors in both)
    CHECK_NEAR(CGAL::to_double((b - a) * (c - a)), (nb - na) * (nc - na), tol);

    // squared area, squared distance, centroid
    CHECK_NEAR(CGAL::to_double(CGAL::squared_area(a, b, c)),
               nat::squared_area(na, nb, nc), tol);
    CHECK_NEAR(CGAL::to_double(CGAL::squared_distance(a, b)),
               nat::squared_distance(na, nb), tol);
    const CP cc = CGAL::centroid(a, b, c);
    const nat::Point3 ncc = nat::centroid(na, nb, nc);
    CHECK_NEAR(CGAL::to_double(cc.x()), ncc.x(), tol);
    CHECK_NEAR(CGAL::to_double(cc.y()), ncc.y(), tol);
    CHECK_NEAR(CGAL::to_double(cc.z()), ncc.z(), tol);
  }
}
