// Tier 3 — local unit tests pinning each internal component against an
// independent oracle (analytic value or finite difference). Written on the CGAL
// backend; the same tests validate the native backend once the seam lands.
#include "test_framework.h"
#include "test_backend.h"

#include "convex_triangulation.h"
#include "dual_hull.h"
#include "polyhedron_utils.h"
#include "vertex_erasure.h"

#include <igl/icosahedron.h>

#include <CGAL/convex_hull_3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

namespace
{
using K = TestPolyhedron::Traits;   // Epick
using P3 = TestPolyhedron::Point_3;

// Build a convex-hull polyhedron from Eigen points, with contiguous v/f ids.
TestPolyhedron build_hull(const Eigen::MatrixXd & V)
{
  auto pts = point_list<K>(V);
  TestPolyhedron poly;
  CGAL::convex_hull_3(pts.begin(), pts.end(), poly);
  int vid = 0;
  for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) v->id() = vid++;
  int fid = 0;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f) f->id() = fid++;
  return poly;
}

// Reproduce ConvexHullSimplification's primal->dual construction for a clean
// (degeneracy-free) input. Returns the dual polyhedron and Chebyshev center.
std::pair<TestPolyhedron, P3> build_dual(const Eigen::MatrixXd & V)
{
  TestPolyhedron primal = build_hull(V);
  const double tol = 1e-20;
  P3 x0 = polyhedron_chebyshev_center(primal, tol);
  auto dpts = dual_points_list(primal, tol, x0);
  Eigen::MatrixXd dV(dpts.size(), 3);
  for(int i = 0; i < (int)dpts.size(); ++i)
    dV.row(i) << CGAL::to_double(dpts[i].x()), CGAL::to_double(dpts[i].y()),
                 CGAL::to_double(dpts[i].z());
  TestPolyhedron dual = build_hull(dV);
  return {dual, x0};
}

double bbox_diag(const Eigen::MatrixXd & V)
{
  return (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
}

// Euler characteristic V-E+F of a closed pure-triangle polyhedron.
int euler(const TestPolyhedron & p)
{
  return int(p.size_of_vertices()) - int(p.size_of_halfedges() / 2) +
         int(p.size_of_facets());
}

// Signed volume of a polygon mesh (fan-triangulated from each polygon's first
// vertex; tetra volumes from the origin).
double polygon_mesh_volume(const Eigen::MatrixXd & pV,
                           const Eigen::VectorXi & pPI,
                           const Eigen::VectorXi & pPC)
{
  double vol = 0.0;
  for(int f = 0; f + 1 < pPC.size(); ++f)
  {
    const int s = pPC(f), n = pPC(f + 1) - pPC(f);
    const Eigen::RowVector3d a = pV.row(pPI(s));
    for(int i = 1; i + 1 < n; ++i)
    {
      const Eigen::RowVector3d b = pV.row(pPI(s + i));
      const Eigen::RowVector3d c = pV.row(pPI(s + i + 1));
      vol += a.dot(b.cross(c)) / 6.0;
    }
  }
  return std::fabs(vol);
}

// Surface area of a polygon mesh via Newell's method per polygon.
double polygon_mesh_area(const Eigen::MatrixXd & pV,
                         const Eigen::VectorXi & pPI,
                         const Eigen::VectorXi & pPC)
{
  double area = 0.0;
  for(int f = 0; f + 1 < pPC.size(); ++f)
  {
    const int s = pPC(f), n = pPC(f + 1) - pPC(f);
    Eigen::RowVector3d nrm = Eigen::RowVector3d::Zero();
    for(int i = 0; i < n; ++i)
    {
      const Eigen::RowVector3d a = pV.row(pPI(s + i));
      const Eigen::RowVector3d b = pV.row(pPI(s + (i + 1) % n));
      nrm += a.cross(b);
    }
    area += 0.5 * nrm.norm();
  }
  return area;
}
}  // namespace

// --- Chebyshev center -----------------------------------------------------

TEST(comp_chebyshev_cube)
{
  TestPolyhedron cube = build_hull(cube_vertices());
  P3 c = polyhedron_chebyshev_center(cube, 1e-20);
  CHECK_NEAR(CGAL::to_double(c.x()), 0.0, 1e-9);
  CHECK_NEAR(CGAL::to_double(c.y()), 0.0, 1e-9);
  CHECK_NEAR(CGAL::to_double(c.z()), 0.0, 1e-9);
}

TEST(comp_chebyshev_tetrahedron)
{
  TestPolyhedron tet = build_hull(regular_tetrahedron_vertices());
  P3 c = polyhedron_chebyshev_center(tet, 1e-20);
  // Centroid of the regular tetrahedron is the origin.
  CHECK_NEAR(CGAL::to_double(c.x()), 0.0, 1e-9);
  CHECK_NEAR(CGAL::to_double(c.y()), 0.0, 1e-9);
  CHECK_NEAR(CGAL::to_double(c.z()), 0.0, 1e-9);
}

// --- Dual polarity round-trip --------------------------------------------

TEST(comp_dual_roundtrip_icosahedron)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);

  auto [dual, x0] = build_dual(V);
  // dual faces -> primal vertices. The hull triangulates the dual's coplanar
  // (pentagonal) faces, so there are more recovered points than icosahedron
  // vertices, but every icosahedron vertex must be among them.
  auto pV = dual_to_primal_points(dual, x0);
  CHECK(pV.rows() >= V.rows());

  const double tol = 1e-9 * bbox_diag(V);
  for(int i = 0; i < V.rows(); ++i)
  {
    double best = 1e30;
    for(int j = 0; j < pV.rows(); ++j)
      best = std::min(best,
        (V.row(i) - Eigen::RowVector3d(CGAL::to_double(pV(j, 0)),
                                       CGAL::to_double(pV(j, 1)),
                                       CGAL::to_double(pV(j, 2)))).norm());
    CHECK_MSG(best < tol, "primal vertex " + std::to_string(i) +
              " not recovered by dual round-trip");
  }
}

// --- primal_change costs vs finite-difference oracle ----------------------

TEST(comp_primal_change_volume_finite_difference)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  ConvexHullSimplification chs(V, F, 100, CostFunction::PRIMAL_VOLUME);

  auto [pV0, pPI0, pPC0] = chs.get_primal_mesh();
  const double vol0 = polygon_mesh_volume(pV0, pPI0, pPC0);

  CHECK(chs.step());  // remove the cheapest vertex
  const double reported = chs.popped_dual_vertex_costs()(chs.popped_dual_vertex_ids()(0));

  auto [pV1, pPI1, pPC1] = chs.get_primal_mesh();
  const double vol1 = polygon_mesh_volume(pV1, pPI1, pPC1);

  CHECK_NEAR(vol1 - vol0, reported, 1e-9 * vol0);
}

TEST(comp_primal_change_area_finite_difference)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  ConvexHullSimplification chs(V, F, 100, CostFunction::PRIMAL_AREA);

  auto [pV0, pPI0, pPC0] = chs.get_primal_mesh();
  const double a0 = polygon_mesh_area(pV0, pPI0, pPC0);

  CHECK(chs.step());
  const double reported = chs.popped_dual_vertex_costs()(chs.popped_dual_vertex_ids()(0));

  auto [pV1, pPI1, pPC1] = chs.get_primal_mesh();
  const double a1 = polygon_mesh_area(pV1, pPI1, pPC1);

  // area cost is the (positive) change in surface area.
  CHECK_NEAR(a1 - a0, reported, 1e-9 * a0);
}

// --- mean width vs analytic ----------------------------------------------

TEST(comp_mean_width_cube)
{
  Eigen::MatrixXd V = cube_vertices();  // [-1,1]^3, edge length 2
  Eigen::MatrixXi F(0, 3);
  ConvexHullSimplification chs(V, F);
  // Mean width of a box a x b x c is (a+b+c)/2; cube edge 2 -> 3.0.
  CHECK_NEAR(chs.mean_width(), 3.0, 1e-9);
}

TEST(comp_mean_width_icosahedron_regression)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  ConvexHullSimplification chs(V, F);
  // Pinned from the working CGAL baseline.
  CHECK_NEAR(chs.mean_width(), 1.832, 1e-3);
}

// --- half-edge Euler operations ------------------------------------------

TEST(comp_euler_ops_octahedron)
{
  Eigen::MatrixXd V(6, 3);
  V << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
  TestPolyhedron oct = build_hull(V);
  CHECK(oct.is_pure_triangle());
  CHECK(oct.is_closed());
  CHECK(oct.size_of_vertices() == 6);
  CHECK(oct.size_of_facets() == 8);
  CHECK(oct.size_of_halfedges() == 24);
  CHECK(euler(oct) == 2);

  // flip_edge: an interior edge flip preserves counts and closedness.
  {
    auto h = oct.halfedges_begin();
    while(h->is_border_edge()) ++h;
    CGAL::Euler::flip_edge(h, oct);
    CHECK(oct.is_closed());
    CHECK(oct.size_of_vertices() == 6);
    CHECK(oct.size_of_facets() == 8);
    CHECK(euler(oct) == 2);
  }

  // erase_center_vertex on a degree-4 vertex: star of 4 tris -> 1 quad facet.
  {
    TestPolyhedron oct2 = build_hull(V);
    auto v = oct2.vertices_begin();
    const int deg = v->degree();
    auto g = oct2.erase_center_vertex(v->halfedge());
    CHECK(oct2.is_closed());
    CHECK(oct2.size_of_vertices() == 5);
    CHECK(oct2.size_of_facets() == 8u - deg + 1);
    CHECK(euler(oct2) == 2);
    CHECK(!oct2.is_pure_triangle());  // the merged facet is a polygon

    // split_facet re-triangulates the merged facet.
    auto a = g;
    auto b = a->next()->next();
    oct2.split_facet(a, b);
    CHECK(oct2.is_closed());
    CHECK(euler(oct2) == 2);
  }
}

// --- convex one-ring triangulation ---------------------------------------

TEST(comp_one_ring_apex_on_negative_side)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  auto [dual, x0] = build_dual(V);

  // The removed apex must NOT be on the negative side of any new face — this is
  // exactly the convexity invariant measure_vertex_erasure asserts (the apex
  // sits outward of the convex fill). Generically strictly positive; coplanar
  // configurations land on the boundary.
  auto apex_not_negative = [](TestPolyhedron & ring, const P3 & apex)
  {
    for(auto f = ring.facets_begin(); f != ring.facets_end(); ++f)
    {
      auto h = f->halfedge();
      CHECK(CGAL::orientation(h->vertex()->point(),
                              h->next()->vertex()->point(),
                              h->next()->next()->vertex()->point(),
                              apex) != CGAL::NEGATIVE);
    }
  };
  for(auto v = dual.vertices_begin(); v != dual.vertices_end(); ++v)
  {
    const P3 apex = v->point();
    const int deg = v->degree();

    TestPolyhedron ring_flips, ring_hull;
    TestPolyhedron::Halfedge_handle h0a, h0b;
    one_ring_triangulation_convex_via_flips(dual, v, ring_flips, h0a);
    one_ring_triangulation_convex_via_convex_hull(dual, v, ring_hull, h0b);

    // Core convexity invariant for BOTH methods.
    apex_not_negative(ring_flips, apex);
    apex_not_negative(ring_hull, apex);

    // The ring is the triangulated link polygon (the deg neighbors, without the
    // apex). The flip method preserves every neighbor (it only flips edges): deg
    // vertices, deg-2 triangles. The convex-hull method may legitimately drop
    // recessed neighbors (interior to the hull of the others), so it has <= deg.
    CHECK(int(ring_flips.size_of_vertices()) == deg);
    CHECK(int(ring_flips.size_of_facets()) == deg - 2);
    CHECK(int(ring_hull.size_of_vertices()) <= deg);
  }
}

// --- clip_ears path replay fidelity --------------------------------------

TEST(comp_clip_ears_replay_matches_measurement)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::icosahedron(V, F);
  auto [dual, x0] = build_dual(V);

  // Measure erasure of one vertex: yields the recorded ear-clipping path plus
  // (independently) the convex one-ring triangulation it corresponds to.
  auto v = dual.vertices_begin();
  const int vid = v->id();

  // The triangle set the measurement's convex one-ring produces.
  TestPolyhedron ring;
  TestPolyhedron::Halfedge_handle h0;
  one_ring_triangulation_convex_via_flips(dual, v, ring, h0);
  std::set<std::array<int, 3>> expected;
  for(auto f = ring.facets_begin(); f != ring.facets_end(); ++f)
  {
    auto h = f->halfedge();
    std::array<int, 3> t = {int(h->vertex()->id()),
                            int(h->next()->vertex()->id()),
                            int(h->next()->next()->vertex()->id())};
    std::sort(t.begin(), t.end());
    expected.insert(t);
  }

  auto [cost, record] = measure_vertex_erasure(dual, v, 100, CostFunction::PRIMAL_VOLUME);

  // Replay the recorded path on the actual dual.
  auto vh = v;
  erase_vertex_and_clip_ears(dual, vh, record);
  CHECK(dual.is_pure_triangle());
  CHECK(dual.is_closed());
  CHECK(euler(dual) == 2);

  // The faces that fill the hole (all faces touching only the old neighbors, not
  // vid) must be exactly the measured convex triangulation.
  std::set<std::array<int, 3>> got;
  for(auto f = dual.facets_begin(); f != dual.facets_end(); ++f)
  {
    auto h = f->halfedge();
    std::array<int, 3> t = {int(h->vertex()->id()),
                            int(h->next()->vertex()->id()),
                            int(h->next()->next()->vertex()->id())};
    if(t[0] == vid || t[1] == vid || t[2] == vid) continue;
    std::sort(t.begin(), t.end());
    if(expected.count(t)) got.insert(t);
  }
  CHECK(got == expected);
}
