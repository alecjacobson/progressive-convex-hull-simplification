// Phase B, brick 2 — validate the native half-edge mesh: build, invariants,
// flip_edge, and the retriangulate_star splice.
#include "test_framework.h"

#include "native_mesh.h"

#include <algorithm>
#include <array>
#include <vector>

namespace
{
using nat::Mesh;
using nat::Point3;

// Octahedron with outward-consistent winding (origin strictly inside every
// face). Any input winding is fixed to outward here.
void make_octahedron(std::vector<Point3> & pts,
                     std::vector<std::array<int, 3>> & faces)
{
  pts = {Point3(1, 0, 0),  Point3(-1, 0, 0), Point3(0, 1, 0),
         Point3(0, -1, 0), Point3(0, 0, 1),  Point3(0, 0, -1)};
  faces = {{0, 2, 4}, {0, 4, 3}, {0, 3, 5}, {0, 5, 2},
           {1, 4, 2}, {1, 3, 4}, {1, 5, 3}, {1, 2, 5}};
  const Point3 o(0, 0, 0);
  for(auto & f : faces)
    if(nat::orient3d(pts[f[0]], pts[f[1]], pts[f[2]], o) == nat::POSITIVE)
      std::swap(f[1], f[2]);   // make origin the negative (inside) side
}

// Consistency check: every face has the origin strictly inside (NEGATIVE side).
void check_outward_closed(Mesh & M, const Point3 & inside)
{
  CHECK(M.is_closed());
  CHECK(M.is_pure_triangle());
  CHECK(M.euler() == 2);
  for(auto f = M.facets_begin(); f != M.facets_end(); ++f)
  {
    auto h = f->halfedge();
    CHECK(nat::orient3d(h->vertex()->point(),
                        h->next()->vertex()->point(),
                        h->next()->next()->vertex()->point(),
                        inside) == nat::NEGATIVE);
  }
}
}  // namespace

TEST(native_mesh_build_octahedron)
{
  std::vector<Point3> pts;
  std::vector<std::array<int, 3>> faces;
  make_octahedron(pts, faces);
  Mesh M; M.build(pts, faces);

  CHECK(M.size_of_vertices() == 6);
  CHECK(M.size_of_facets() == 8);
  CHECK(M.size_of_halfedges() == 24);   // 12 edges * 2
  check_outward_closed(M, Point3(0, 0, 0));

  // Every octahedron vertex has degree 4.
  for(auto v = M.vertices_begin(); v != M.vertices_end(); ++v)
    CHECK(v->degree() == 4);

  // Opposite pairing and next/prev consistency.
  for(auto e = M.edges_begin(); e != M.edges_end(); ++e)
  {
    CHECK(e->opposite()->opposite() == e);
    CHECK(e->next()->prev() == e);
    CHECK(e->next()->next()->next() == e);   // triangle
  }
}

TEST(native_mesh_open_triangle_has_border)
{
  std::vector<Point3> pts = {Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 1, 0)};
  std::vector<std::array<int, 3>> faces = {{0, 1, 2}};
  Mesh M; M.build(pts, faces);
  CHECK(M.size_of_vertices() == 3);
  CHECK(M.size_of_facets() == 1);
  CHECK(M.size_of_halfedges() == 6);   // 3 interior + 3 border
  CHECK(!M.is_closed());
  CHECK(M.euler() == 1);               // a disk
  // The three border half-edges form a loop.
  int border = 0;
  for(auto h = M.halfedges_begin(); h != M.halfedges_end(); ++h)
    if(h->is_border()) { ++border; CHECK(h->next()->is_border()); }
  CHECK(border == 3);
}

TEST(native_mesh_flip_edge)
{
  std::vector<Point3> pts;
  std::vector<std::array<int, 3>> faces;
  make_octahedron(pts, faces);
  Mesh M; M.build(pts, faces);

  // Flip one interior edge; structure and closedness must be preserved.
  auto h = M.halfedges_begin();
  while(h->is_border_edge()) ++h;
  M.flip_edge(h);
  CHECK(M.is_closed());
  CHECK(M.is_pure_triangle());
  CHECK(M.euler() == 2);
  CHECK(M.size_of_vertices() == 6);
  CHECK(M.size_of_facets() == 8);
  for(auto e = M.edges_begin(); e != M.edges_end(); ++e)
  {
    CHECK(e->next()->next()->next() == e);
    CHECK(e->opposite()->opposite() == e);
  }
}

TEST(native_mesh_retriangulate_star)
{
  std::vector<Point3> pts;
  std::vector<std::array<int, 3>> faces;
  make_octahedron(pts, faces);
  Mesh M; M.build(pts, faces);

  auto v0 = M.vertices_begin();          // vertex (1,0,0), degree 4
  const Point3 apex = v0->point();
  auto ring = M.ring_vertices(v0);
  const int k = (int)ring.size();
  CHECK(k == 4);
  std::vector<int> ring_idx;
  for(auto r : ring) ring_idx.push_back(r.index());
  auto is_ring = [&](int i){ return std::find(ring_idx.begin(), ring_idx.end(), i) != ring_idx.end(); };

  // Fan triangulation of the ring polygon: (0, i, i+1).
  std::vector<std::array<int, 3>> tris;
  for(int i = 1; i + 1 < k; ++i) tris.push_back({0, i, i + 1});

  M.retriangulate_star(v0, tris);

  CHECK(M.size_of_vertices() == 5);
  CHECK(M.size_of_facets() == 6);        // 8 - 4 + 2
  CHECK(M.is_closed());
  CHECK(M.is_pure_triangle());
  CHECK(M.euler() == 2);

  // The removed apex must not be on the negative side of any NEW fill face
  // (same convexity invariant the algorithm relies on). New faces are exactly
  // those whose three vertices are all ring vertices; far-side original faces
  // legitimately keep the apex on their interior/negative side.
  int fill_faces = 0;
  for(auto f = M.facets_begin(); f != M.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const int a = h->vertex().index();
    const int b = h->next()->vertex().index();
    const int c = h->next()->next()->vertex().index();
    if(!(is_ring(a) && is_ring(b) && is_ring(c))) continue;
    ++fill_faces;
    CHECK(nat::orient3d(h->vertex()->point(),
                        h->next()->vertex()->point(),
                        h->next()->next()->vertex()->point(),
                        apex) != nat::NEGATIVE);
  }
  CHECK(fill_faces == k - 2);   // a triangulated k-gon

  // garbage_collect keeps it valid.
  M.garbage_collect();
  CHECK(M.size_of_vertices() == 5);
  CHECK(M.size_of_facets() == 6);
  CHECK(M.is_closed());
  CHECK(M.euler() == 2);
}
