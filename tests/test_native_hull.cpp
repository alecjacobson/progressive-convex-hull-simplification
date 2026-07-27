// Phase B, brick 3 — validate nat::convex_hull_3 (qhull) against CGAL on shared
// point sets: same hull-vertex count, containment of all inputs, closed/euler,
// and matching volume.
#include "test_framework.h"

#include "native_hull.h"
#include "test_backend.h"        // TestPolyhedron, cube/tetra helpers

#include "polyhedron_utils.h"

#include <igl/icosahedron.h>

#include <CGAL/convex_hull_3.h>

#include <array>
#include <cmath>
#include <vector>

namespace
{
using nat::Mesh;
using nat::Point3;

std::vector<Point3> to_native(const Eigen::MatrixXd & V)
{
  std::vector<Point3> pts(V.rows());
  for(int i = 0; i < V.rows(); ++i) pts[i] = Point3(V(i, 0), V(i, 1), V(i, 2));
  return pts;
}

// Signed volume of a closed native triangle mesh (tetrahedra from origin).
double mesh_volume(Mesh & M)
{
  double vol = 0.0;
  for(auto f = M.facets_begin(); f != M.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const Point3 a = h->vertex()->point();
    const Point3 b = h->next()->vertex()->point();
    const Point3 c = h->next()->next()->vertex()->point();
    vol += (a - nat::ORIGIN) * nat::cross_product(b - nat::ORIGIN, c - nat::ORIGIN) / 6.0;
  }
  return std::fabs(vol);
}

double cgal_hull_volume(const Eigen::MatrixXd & V, int & nverts)
{
  using K = TestPolyhedron::Traits;
  auto pts = point_list<K>(V);
  TestPolyhedron poly;
  CGAL::convex_hull_3(pts.begin(), pts.end(), poly);
  nverts = (int)poly.size_of_vertices();
  double vol = 0.0;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    auto a = h->vertex()->point();
    auto b = h->next()->vertex()->point();
    auto c = h->next()->next()->vertex()->point();
    vol += CGAL::to_double((a - CGAL::ORIGIN) *
             CGAL::cross_product(b - CGAL::ORIGIN, c - CGAL::ORIGIN)) / 6.0;
  }
  return std::fabs(vol);
}

// Assert every input point is inside (within tol) the native hull M, using a
// tolerant signed-distance test (points that ARE hull vertices sit exactly on
// their incident faces, so exact orient3d sign is too brittle here).
void check_contains(Mesh & M, const std::vector<Point3> & pts, double scale)
{
  const double tol = 1e-7 * scale;
  for(const auto & p : pts)
    for(auto f = M.facets_begin(); f != M.facets_end(); ++f)
    {
      auto h = f->halfedge();
      const Point3 a = h->vertex()->point();
      const Point3 b = h->next()->vertex()->point();
      const Point3 c = h->next()->next()->vertex()->point();
      const nat::Vector3 n = nat::cross_product(b - a, c - a);  // outward
      const double len = std::sqrt(n * n);
      if(len < 1e-20) continue;
      const double dist = (n * (p - a)) / len;  // >0 means outside
      CHECK_MSG(dist <= tol, "input point outside hull by " + std::to_string(dist));
    }
}
}  // namespace

TEST(native_hull_cube_with_interior)
{
  // 8 corners + interior points; hull must be the 8-corner cube.
  Eigen::MatrixXd V(11, 3);
  V << -1,-1,-1,  1,-1,-1,  -1,1,-1,  1,1,-1,
       -1,-1, 1,  1,-1, 1,  -1,1, 1,  1,1, 1,
        0, 0, 0,  0.5,0,0,  0,-0.5,0.3;
  Mesh M = nat::convex_hull_3(to_native(V));
  CHECK(M.size_of_vertices() == 8);
  CHECK(M.size_of_facets() == 12);   // triangulated cube
  CHECK(M.is_closed());
  CHECK(M.is_pure_triangle());
  CHECK(M.euler() == 2);
  CHECK_NEAR(mesh_volume(M), 8.0, 1e-9);
}

TEST(native_hull_icosahedron_matches_cgal)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  auto pts = to_native(V);
  Mesh M = nat::convex_hull_3(pts);

  CHECK(M.is_closed());
  CHECK(M.is_pure_triangle());
  CHECK(M.euler() == 2);
  CHECK(M.size_of_vertices() == 12);
  CHECK(M.size_of_facets() == 20);

  int cgal_nv = 0;
  const double cgal_vol = cgal_hull_volume(V, cgal_nv);
  CHECK((int)M.size_of_vertices() == cgal_nv);
  CHECK_NEAR(mesh_volume(M), cgal_vol, 1e-9 * cgal_vol);

  const double scale = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  check_contains(M, pts, scale);
}

TEST(native_hull_random_cloud_matches_cgal)
{
  // Deterministic pseudo-random cloud (LCG) — no RNG headers.
  const int n = 300;
  Eigen::MatrixXd V(n, 3);
  unsigned long s = 123456789UL;
  auto rnd = [&]() { s = s * 6364136223846793005ULL + 1442695040888963407ULL;
                     return ((s >> 33) & 0xFFFFFF) / double(0xFFFFFF) * 2.0 - 1.0; };
  for(int i = 0; i < n; ++i) V.row(i) << rnd(), rnd(), rnd();

  auto pts = to_native(V);
  Mesh M = nat::convex_hull_3(pts);
  CHECK(M.is_closed());
  CHECK(M.is_pure_triangle());
  CHECK(M.euler() == 2);

  int cgal_nv = 0;
  const double cgal_vol = cgal_hull_volume(V, cgal_nv);
  CHECK((int)M.size_of_vertices() == cgal_nv);
  CHECK_NEAR(mesh_volume(M), cgal_vol, 1e-9 * cgal_vol);
  check_contains(M, pts, 2.0);
}
