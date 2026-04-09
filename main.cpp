#include "chebyshev_center.h"
#include "primal_volume_subtended.hpp"

#include <igl/read_triangle_mesh.h>
#include <igl/min_heap.h>
#include <igl/icosahedron.h>
#include <igl/matlab_format.h>
#include <igl/writeDMAT.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/get_seconds.h>
#include <igl/copyleft/cgal/assign.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/join_coplanar_neighboring_facets.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_copy_3.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>
#include <Eigen/Core>
#include <iostream>
#include <cstdio>

std::vector<int> num_flips_list;
std::vector<int> local_size_list;
std::vector<double> flip_times_list;

int MAX_DEGREE_FOR_FLIPS = 100;

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


template <typename Kernel>
Eigen::Matrix<typename Kernel::FT,Eigen::Dynamic,3,Eigen::RowMajor> point_matrix(
    const std::vector<typename Kernel::Point_3> & points)
{
  Eigen::Matrix<typename Kernel::FT,Eigen::Dynamic,3,Eigen::RowMajor> V(points.size(),3);
  for(int i = 0;i<V.rows();i++)
  {
    V(i,0) = points[i].x();
    V(i,1) = points[i].y();
    V(i,2) = points[i].z();
  }
  return V;
}

// Given a triangulated CGAL polyhedron, compute the Chebyshev center of the
// halfspaces corresponding to faces with squared area > primal_squared_area_tol.
// Converts to double to call igl::chebyshev_center, returns result in Kernel.
template <class Polyhedron>
typename Polyhedron::Traits::Point_3 chebyshev_center(
  const Polyhedron & poly,
  const double primal_squared_area_tol)
{
  using Kernel = typename Polyhedron::Traits;
  using Point = typename Kernel::Point_3;

  // Two-pass: count then fill
  int nf = 0;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    if(CGAL::to_double(CGAL::squared_area(
        h->vertex()->point(),
        h->next()->vertex()->point(),
        h->next()->next()->vertex()->point())) > primal_squared_area_tol) { nf++; }
  }
  Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor> P(nf, 4);
  int i = 0;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const auto & A = h->vertex()->point();
    const auto & B = h->next()->vertex()->point();
    const auto & C = h->next()->next()->vertex()->point();
    if(CGAL::to_double(CGAL::squared_area(A,B,C)) <= primal_squared_area_tol) { continue; }
    const Eigen::Vector3d a(CGAL::to_double(A.x()),CGAL::to_double(A.y()),CGAL::to_double(A.z()));
    const Eigen::Vector3d b(CGAL::to_double(B.x()),CGAL::to_double(B.y()),CGAL::to_double(B.z()));
    const Eigen::Vector3d c(CGAL::to_double(C.x()),CGAL::to_double(C.y()),CGAL::to_double(C.z()));
    const Eigen::Vector3d n = (b-a).cross(c-a);
    const Eigen::Vector3d bc = (a+b+c)/3.0;
    P.row(i++) << n(0), n(1), n(2), -n.dot(bc);
  }
  Eigen::Matrix<double,3,1> x0d;
  igl::chebyshev_center(P, x0d);
  return Point(x0d(0), x0d(1), x0d(2));
}

// Given a triangulated CGAL polyhedron and interior point x0, return one dual
// point per face (skipping faces with squared area <= primal_squared_area_tol).
// Duality map: given halfspace [n, b] (b = -n·bc), d = -n / (n·x0 + b).
// All arithmetic stays in Kernel.
template <class Polyhedron>
std::vector<typename Polyhedron::Traits::Point_3> dual_points_list(
  const Polyhedron & poly,
  const double primal_squared_area_tol,
  const typename Polyhedron::Traits::Point_3 & x0)
{
  using Kernel = typename Polyhedron::Traits;
  using FT = typename Kernel::FT;
  using Point = typename Kernel::Point_3;
  using Vector = typename Kernel::Vector_3;

  std::vector<Point> dpts;
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const auto & A = h->vertex()->point();
    const auto & B = h->next()->vertex()->point();
    const auto & C = h->next()->next()->vertex()->point();
    if(CGAL::to_double(CGAL::squared_area(A,B,C)) <= primal_squared_area_tol) { continue; }
    const Vector n = CGAL::cross_product(B - A, C - A);
    const Vector bc = ((A - CGAL::ORIGIN) + (B - CGAL::ORIGIN) + (C - CGAL::ORIGIN)) / FT(3);
    const FT b = -(n * bc);
    const FT denom = n * (x0 - CGAL::ORIGIN) + b;
    dpts.push_back(CGAL::ORIGIN + (-(n / denom)));
  }
  return dpts;
}

template <class Polyhedron>
void print_in_matlab_format( const Polyhedron & poly,const std::string prefix = "")
{
  Eigen::Matrix<typename Polyhedron::Traits::FT,Eigen::Dynamic,3,Eigen::RowMajor> V;
  Eigen::VectorXi PI,PC;
  igl::copyleft::cgal::polyhedron_to_mesh(poly,V,PI,PC);
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> double_V;
  igl::copyleft::cgal::assign(V,double_V);
  std::cout<<igl::matlab_format(double_V,prefix+"V")<<std::endl;
  std::cout<<igl::matlab_format_index(PI.transpose().eval(),prefix+"PI")<<std::endl;
  std::cout<<igl::matlab_format(      PC.transpose().eval(),      prefix+"PC")<<std::endl;
}

template <class Polyhedron>
void print_VF(const Polyhedron & poly, const std::string prefix = "")
{
  int max_id = -1;
  printf("print_VF\n");
  for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
  {
    printf("v->id() = %d\n",v->id());
    if(int(v->id()) > max_id) { max_id = int(v->id()); }
  }
  printf("-----------------\n");
  printf("max_id = %d\n",max_id);
  Eigen::Matrix<typename Polyhedron::Traits::FT,Eigen::Dynamic,3,Eigen::RowMajor> V(max_id+1,3);
  V.setConstant(std::numeric_limits<typename Polyhedron::Traits::FT>::quiet_NaN());
  for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
  {
    V.row(v->id()) << v->point().x(), v->point().y(), v->point().z();
  }
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
    assert(k == poly.size_of_facets());
  }
  std::cout<<igl::matlab_format(V,prefix+"V")<<std::endl;
  std::cout<<igl::matlab_format_index(F,prefix+"F")<<std::endl;
}

template <class Polyhedron>
void print_faces(const Polyhedron & poly)
{
  for(auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
  {
    auto h = f->halfedge();
    printf("  (%d %d %d)\n",
      h->vertex()->id(),
      h->next()->vertex()->id(),
      h->next()->next()->vertex()->id());
  }
}

template <class Polyhedron>
Polyhedron extract_copy_of_one_ring(
  const Polyhedron& poly,
  typename Polyhedron::Vertex_const_handle v)
{

  using Point = typename Polyhedron::Point_3;
  Polyhedron out;

  auto h0 = v->halfedge();
  if(h0 == nullptr) return out;

  std::vector<Point> points;
  std::vector<std::size_t> vertex_ids;

  points.push_back(v->point());
  vertex_ids.push_back(v->id());

  // collect neighboring vertices in cyclic order around v
  auto h = h0;
  do
  {
    points.push_back(h->opposite()->vertex()->point());
    vertex_ids.push_back(h->opposite()->vertex()->id());
    h = h->next()->opposite();
  } while(h != h0);

  struct Builder : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
  {
    const std::vector<Point>& points;
    Builder(const std::vector<Point>& points) : points(points) {}

    void operator()(typename Polyhedron::HalfedgeDS& hds)
    {
      CGAL::Polyhedron_incremental_builder_3<typename Polyhedron::HalfedgeDS> B(hds, true);
      const int n = static_cast<int>(points.size()) - 1;

      B.begin_surface(points.size(), n);

      for(const auto& p : points)
      {
        B.add_vertex(p);
      }

      for(int i = 1; i <= n; ++i)
      {
        const int j = (i < n ? i+1 : 1);
        B.begin_facet();
        B.add_vertex_to_facet(0);
        B.add_vertex_to_facet(j);
        B.add_vertex_to_facet(i);
        B.end_facet();
      }

      B.end_surface();
    }
  };

  Builder builder(points);
  out.delegate(builder);

  // Copy vertex IDs in insertion order
  {
    auto vv = out.vertices_begin();
    for(std::size_t i = 0; i < vertex_ids.size(); ++i, ++vv)
    {
      vv->id() = vertex_ids[i];
    }
  }

  return out;
}

template <class Polyhedron>
typename Polyhedron::Traits::Oriented_side 
edge_side(
  typename Polyhedron::Halfedge_const_handle e)
{
  using Point = typename Polyhedron::Point_3;
  using Halfedge_handle = typename Polyhedron::Halfedge_const_handle;

  Halfedge_handle h1 = e;
  Halfedge_handle h2 = e->opposite();

  const Point& a1 = h1->vertex()->point();
  const Point& b1 = h1->next()->vertex()->point();
  const Point& c1 = h1->next()->next()->vertex()->point();
  // vertex of f2 opposite the shared edge h1
  const Point& p2 = h2->next()->vertex()->point();

  // Don't do this! This counts as a construction!
  //Plane P1(a1,b1,c1);
  //auto side = P1.oriented_side(p2);
  
  const auto ori = CGAL::orientation(a1,b1,c1,p2);
  switch(ori)
  {
    case CGAL::POSITIVE: return CGAL::ON_POSITIVE_SIDE;
    case CGAL::NEGATIVE: return CGAL::ON_NEGATIVE_SIDE;
    case CGAL::COPLANAR: return CGAL::ON_ORIENTED_BOUNDARY;
  }
}

template <typename Polyhedron>
void confirm_all_edges_are_convex(const Polyhedron& poly)
{
  using Halfedge_handle = typename Polyhedron::Halfedge_const_handle;
  using Point           = typename Polyhedron::Point_3;

  // Map all vertex handles to ids
  std::map<typename Polyhedron::Vertex_const_handle, int> vertex_id_map;
  {
    int vid = 0;
    for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    {
      vertex_id_map[v] = vid++;
    }
  }

  for(auto e = poly.halfedges_begin(); e != poly.halfedges_end(); ++e)
  {
    // process each undirected edge once
    if(e > e->opposite()) continue;

    if(e->is_border_edge()) continue;

    auto side = edge_side<Polyhedron>(e);

    //{
    //  auto ee = e;
    //  auto side = edge_side<Polyhedron>(ee);
    //  Halfedge_handle h1 = ee;
    //  Halfedge_handle h2 = ee->opposite();
    //  printf("(%d,%d) → (%d,%d,%d) + %d → %s\n",
    //    vertex_id_map[h1->vertex()], vertex_id_map[h2->vertex()],
    //    vertex_id_map[h1->vertex()], vertex_id_map[h1->next()->vertex()], vertex_id_map[h1->next()->next()->vertex()],
    //    vertex_id_map[h2->next()->vertex()],
    //    (side == CGAL::ON_NEGATIVE_SIDE ? "convex" : (side == CGAL::ON_POSITIVE_SIDE ? "non-convex" : "coplanar"))
    //  );
    //}
    //{
    //  auto ee = e->opposite();
    //  auto side = edge_side<Polyhedron>(ee);
    //  Halfedge_handle h1 = ee;
    //  Halfedge_handle h2 = ee->opposite();
    //  printf("(%d,%d) → (%d,%d,%d) + %d → %s\n",
    //    vertex_id_map[h1->vertex()], vertex_id_map[h2->vertex()],
    //    vertex_id_map[h1->vertex()], vertex_id_map[h1->next()->vertex()], vertex_id_map[h1->next()->next()->vertex()],
    //    vertex_id_map[h2->next()->vertex()],
    //    (side == CGAL::ON_NEGATIVE_SIDE ? "convex" : (side == CGAL::ON_POSITIVE_SIDE ? "non-convex" : "coplanar"))
    //  );
    //}

    if(side == CGAL::ON_POSITIVE_SIDE)
    {
      throw std::runtime_error("Found a non-convex edge.");
    }
  }
}


template <class Polyhedron>
int flip_until_all_interior_edges_are_convex(
  Polyhedron& poly,
  const int max_flips = 1000)
{
#if false

#warning This might even be O(N³). I think the worst case number of flips is \
O(N²), but this implementation restarts the search for a non-convex edge \
after every flip, which could lead to O(N³) behavior. Does the inner for \
loop really need the `break` after the flip, or could it continue to \
halfedges_end() even if it misses some newly non-convex edges knowing \
that it will get them on the next pass?
  int num_flips = 0;
  bool changed = true;
  while(changed)
  {
     changed = false;
     for(auto e = poly.halfedges_begin();
         e != poly.halfedges_end(); ++e)
     {
       // only consider each undirected edge once
       if(e > e->opposite()) continue;
       // only internal edges can be flipped
       if(e->is_border_edge()) continue;
       auto side = edge_side<Polyhedron>(e);
       if(side == CGAL::ON_POSITIVE_SIDE)
       {
         CGAL::Euler::flip_edge(e, poly);
         changed = true;
         num_flips++;
         if(num_flips > max_flips)
         {
           throw std::runtime_error("Too many flips, something is wrong.");
         }
         // Sad O(|E|^2) algorithm. Hoping |E| is small!
         // I claim it's ok to not break. After CGAL::Euler::flip_edge, e is
         // still valid and somewhere in the list of half-edges, so this just
         // continues until halfedges_end() and then starts over at
         // halfedges_begin(), which should be fine. We still potentially
         // recheck everything a full additional time.
         //break;
       }
     }
  }
  return num_flips;
#else
  int num_flips = 0;
  // Push all internal edges onto stack
  std::vector<typename Polyhedron::Halfedge_handle> stack;
  stack.reserve(poly.size_of_halfedges()/2);
  std::function<void(typename Polyhedron::Halfedge_handle)> mark_and_push_edge;
  mark_and_push_edge = [&stack,&mark_and_push_edge](auto e)
  {
    if(e->opposite() < e)
    {
      return mark_and_push_edge(e->opposite());
    }
    if(e->is_border_edge()) return;
    e->id() = false;
    stack.push_back(e);
  };

  for(auto e = poly.halfedges_begin(); e != poly.halfedges_end(); ++e)
  {
    if(e > e->opposite()) continue;
    mark_and_push_edge(e);
  }
  while(!stack.empty())
  {
    auto e = stack.back();
    stack.pop_back();
    // Was marked as convex in the meantime.
    if(e->id()) { continue; }
    auto side = edge_side<Polyhedron>(e);
    if(side == CGAL::ON_POSITIVE_SIDE)
    {
      CGAL::Euler::flip_edge(e, poly);
      num_flips++;
      if(num_flips > max_flips)
      {
        throw std::runtime_error("Too many flips, something is wrong.");
      }
      mark_and_push_edge(e->next());
      mark_and_push_edge(e->next()->next());
      mark_and_push_edge(e->opposite()->next());
      mark_and_push_edge(e->opposite()->next()->next());
    }
  }


#ifndef NDEBUG
     for(auto e = poly.halfedges_begin();
         e != poly.halfedges_end(); ++e)
     {
       // only consider each undirected edge once
       if(e > e->opposite()) continue;
       // only internal edges can be flipped
       if(e->is_border_edge()) continue;
       auto side = edge_side<Polyhedron>(e);
       if(side == CGAL::ON_POSITIVE_SIDE)

        {
          throw std::runtime_error("Found a non-convex edge after flipping, something is wrong.");
        }
     }
#endif

  return num_flips;
#endif
}

template <class Polyhedron>
std::pair< size_t, std::vector<int> >
clip_ears(Polyhedron & one_ring_copy, typename Polyhedron::Halfedge_handle h0)
{
  const auto walk = [](const auto h)->auto
  {
    assert(h->opposite()->is_border());
    return h->opposite()->prev()->opposite();
  };

  // Sanity check.
#ifndef NDEBUG
  {
    int count = 0;
    auto h = h0;
    while(true)
    {
      h = walk(h);
      count++;
      assert(count < 1000);
      if(h == h0){break;}
    }
    assert(count == one_ring_copy.size_of_vertices());
  }
#endif

  int walks_since_last_ear = 0;
  auto h = h0;
  auto start_vertex_id = h->vertex()->id();
  std::vector<int> path;
  path.reserve(one_ring_copy.size_of_facets()-1);
  while(one_ring_copy.size_of_vertices() > 3)
  {
    auto w = walk(h);
    auto n = h->next();
    //printf("h: %d→%d\n",h->opposite()->vertex()->id(),h->vertex()->id());
    if(w == n)
    {
      // one_ring_copy is triangulated so this is dumb to test.
      assert(h->next()->next()->next() == h && "h is on a triangle.");
      // Ear
      //    n₁
      // o--←--o
      //  h   /n₀==w
      //   ↘ ↗
      //    o
      path.push_back(walks_since_last_ear);
      walks_since_last_ear = 0;
      //printf("before clipping:\n");
      //print_faces(one_ring_copy);
      //printf("  clipping ear: (%d %d %d)\n",
      //    h->vertex()->id(),
      //    h->next()->vertex()->id(),
      //    h->next()->next()->vertex()->id());
      auto q = h->next()->next()->opposite();
      one_ring_copy.erase_facet(h);
      assert(q->opposite()->is_border());
      //printf("q: %d→%d\n",q->opposite()->vertex()->id(),q->vertex()->id());
      //printf("after clipping:\n");
      //print_faces(one_ring_copy);
      h = q;
    }else
    {
      // not Ear
      // o--←--o
      //  h  n/  …
      //   ↘ ↗
      //    o-→w-o
      h = w;
      walks_since_last_ear++;
    }
  }
  return {start_vertex_id, path};
}

template <class Polyhedron>
void erase_vertex_and_clip_ears(
  Polyhedron & dual,
  typename Polyhedron::Vertex_handle & v,
  const size_t start_vertex_id,
  const std::vector<int> & path)
{
  auto h0 = dual.erase_center_vertex(v->halfedge());
  while(h0->vertex()->id() != start_vertex_id)
  {
    h0 = h0->next();
  }
  // Now h0 matches h0 on one_ring_copy used to make path.
  const auto walk = [](const auto h)->auto
  {
    return h->next();
  };

  //{
  //  int count = 0;
  //  auto h = h0;
  //  while(true)
  //  {
  //    h = walk(h);
  //    count++;
  //    assert(count < 1000);
  //    if(h == h0){break;}
  //  }
  //  printf("count = %d\n",count);
  //}

  {
    auto h = h0;
    for(int i = 0;i<path.size();i++)
    {
      const int steps = path[i];
      for(int j = 0;j<steps;j++)
      {
        h = walk(h);
      }
      // glue the ear back on
      h = dual.split_facet(h->prev(),h->next());
      //{
      //  auto g = h->opposite();
      //  printf("glue (%d %d %d)\n",
      //      g->vertex()->id(),
      //      g->next()->vertex()->id(),
      //      g->next()->next()->vertex()->id());
      //}
    }
    //{
    //  auto g = h;
    //  printf("left (%d %d %d)\n",
    //      g->vertex()->id(),
    //      g->next()->vertex()->id(),
    //      g->next()->next()->vertex()->id());
    //}
  }
}

struct Record
{
  size_t start_vertex_id;
  std::vector<int> path;
};

// Output by reference so that h0 is valid reference on one_ring_copy.
// Using return value for {one_ring_copy,h0} didn't work.
template <class Polyhedron>
void one_ring_triangulation(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy,
  typename Polyhedron::Halfedge_handle & h0)
{
#if false
  //// print ids of neighbors
  //{
  //  auto h = v->halfedge();
  //  printf("  v id = %d\n", v->id());
  //  printf("    neighbors: ");
  //  do
  //  {
  //    printf("%d ", h->opposite()->vertex()->id());
  //    h = h->next()->opposite();
  //  } while(h != v->halfedge());
  //  printf("\n");
  //}

  // Always in [v, one-ring in order]
  one_ring_copy = extract_copy_of_one_ring(dual, v);
  //// see that ids got copied by printing them
  //{
  //  auto v0 = one_ring_copy.vertices_begin();
  //  printf("  v0 id = %d\n", v0->id());
  //  printf("    neighbors: ");
  //  auto h = v0->halfedge();
  //  do
  //  {
  //    printf("%d ", h->opposite()->vertex()->id());
  //    h = h->next()->opposite();
  //  } while(h != v0->halfedge());
  //  printf("\n");
  //}

  //print_in_matlab_format(dual,"d");
  //print_in_matlab_format(one_ring_copy,"o");

  //printf("confirm_all_edges_are_convex(one_ring_copy)\n");
#ifndef NDEBUG
  confirm_all_edges_are_convex(one_ring_copy);
#endif

  auto v0 = one_ring_copy.vertices_begin();
  // Copy point 
  using Point = typename Polyhedron::Point_3;
  const Point p0 = v0->point();
  // half edge pointing at v0
  auto g = v0->halfedge();
  h0 = one_ring_copy.erase_center_vertex(g);
  //{
  //  printf("  after erasing center vertex, one_ring_copy's facet has vertices: ");
  //  auto h = h0;
  //  while(true)
  //  {
  //    printf("%d ",h->vertex()->id());
  //    h = h->next();
  //    if(h == h0) { break; }
  //  }
  //  printf("\n");
  //}
  {
    auto h = h0;
    // triangulate fan
    while(h->next()->next()->next() != h)
    {
      CGAL::Euler::split_face(h, h->next()->next(), one_ring_copy);
    }
  }
#else
  using Point = typename Polyhedron::Point_3;
  {
    const int nv = v->degree();
     one_ring_copy.reserve(nv, 2*nv + 2*(nv-3), nv-2);
     const auto h_start = v->halfedge()->opposite()->prev();
     auto h = h_start;
     int i = 0;
     std::vector<typename Polyhedron::Vertex_const_handle> prev_verts(2);
     typename Polyhedron::Halfedge_handle q,H;
     do
     {
       //verts[2] = pointer to h->opposite()->vertex()->point();
       const auto vert = h->opposite()->vertex();
       if(i < 2)
       {
         prev_verts[i] = vert;
       }
       else if(i == 2)
       {
         q = one_ring_copy.make_triangle(
           prev_verts[0]->point(),
           prev_verts[1]->point(),
                    vert->point());
         // copy the ids
         q->vertex()->id() = prev_verts[0]->id();
         q->next()->vertex()->id() = prev_verts[1]->id();
         q->next()->next()->vertex()->id() = vert->id();
         // H is always the border halfedge that points to 0
         H = q->next()->opposite();
         assert(H->vertex() == q->vertex());
         assert(H->is_border());
       }else if(i > 2)
       {
         // q verts to 0th vertex
         //
         //   0------3
         //  ↑| ↖\   ↑
         // H||  qG  s
         //  |↓   \↘ |
         //   1----→2
         auto G = q->opposite();
         assert(G->is_border());
         auto s = one_ring_copy.add_vertex_and_facet_to_border(H,G);
         assert(s->opposite()->is_border());
         s->vertex()->point() = vert->point();
         s->vertex()->id() = vert->id();
         q = s->next();
       }
       i++;
       //h = h->next()->opposite();
       h = h->opposite()->prev();
     } while(h != h_start);
    h0 = q;
  }
#endif
}

template <class Polyhedron>
void one_ring_triangulation_convex_via_flips(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy,
  typename Polyhedron::Halfedge_handle & h0)
{
  one_ring_triangulation(dual,v,one_ring_copy,h0);
  //printf(" %d,  %d,  %d\n",one_ring_copy.size_of_vertices(),one_ring_copy.size_of_halfedges(),one_ring_copy.size_of_facets());
  //IGL_TICTOC_LAMBDA;
  //tictoc();
  const int num_flips = flip_until_all_interior_edges_are_convex(one_ring_copy,1000);
  //num_flips_list.push_back(num_flips);
  //local_size_list.push_back(one_ring_copy.size_of_vertices());
  //flip_times_list.push_back(tictoc());
  assert(num_flips < 1000);
  //printf("confirm_all_edges_are_convex(one_ring_copy) after %d flips\n", num_flips);
#ifndef NDEBUG
  confirm_all_edges_are_convex(one_ring_copy);
#endif
}

template <class Polyhedron>
void one_ring_triangulation_convex_via_convex_hull(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  Polyhedron & one_ring_copy,
  typename Polyhedron::Halfedge_handle & h0)
{
  auto p = v->point();

  using Point = typename Polyhedron::Point_3;
  std::map<Point,int> point_to_id;
  std::vector<Point> points;
  points.reserve(v->degree());
  {
    auto h = v->halfedge();
    do
    {
      auto point = h->opposite()->vertex()->point();
      auto id = h->opposite()->vertex()->id();
      point_to_id[point] = id;
      points.push_back(point);
      h = h->next()->opposite();
    } while(h != v->halfedge());
  }
  CGAL::convex_hull_3(points.begin(),points.end(),one_ring_copy);
  for(auto v = one_ring_copy.vertices_begin(); v != one_ring_copy.vertices_end(); ++v)
  {
    auto point = v->point();
    auto id = point_to_id[point];
    v->id() = id;
  }


  while(true)
  {
    bool found_any = false;
    for(auto f = one_ring_copy.facets_begin(); f != one_ring_copy.facets_end(); )
    {
      auto h = f->halfedge();
      f++;
      const auto & a = h->vertex()->point();
      const auto & b = h->next()->vertex()->point();
      const auto & c = h->next()->next()->vertex()->point();
      const auto ori = CGAL::orientation(a,b,c,p);
      if(ori == CGAL::NEGATIVE)
      {
        found_any = true;
        //printf("removing (%d %d %d) because it's oriented negatively with respect to p\n",
        //    h->vertex()->id(),
        //    h->next()->vertex()->id(),
        //    h->next()->next()->vertex()->id());
        one_ring_copy.erase_facet(h);
      }
    }
    if(!found_any){ break; }
  }

  for(auto h = one_ring_copy.halfedges_begin(); h != one_ring_copy.halfedges_end(); ++h)
  {
    if(h->opposite()->is_border() && h->opposite()->vertex()->id() == v->halfedge()->opposite()->vertex()->id())
    {
      h0 = h;
      break;
    }
  }
  assert(one_ring_copy.size_of_vertices() == v->degree());
  //printf("one_ring_copy:\n");
  //print_faces(one_ring_copy);
  //printf("h0: %d→%d\n",h0->opposite()->vertex()->id(),h0->vertex()->id());
  //printf("one_ring_copy:\n");
  //print_faces(one_ring_copy);
  //printf("h0: %d→%d\n",h0->opposite()->vertex()->id(),h0->vertex()->id());
  //printf("\n");
}


template <class Polyhedron>
std::pair<
  typename Polyhedron::Traits::FT,
  Record>
measure_vertex_erasure(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v)
{
  //printf("measure_vertex_erasure(%d)\n", v->id());
  using Scalar = typename Polyhedron::Traits::FT;
#ifndef NDEBUG
  {
    for(auto u = dual.vertices_begin(); u != dual.vertices_end(); ++u)
    {
      if(u->id() == v->id())
      {
        assert(u == v && "vertex ids should be unique");
        break;
      }
    }
  }
#endif

  auto p = v->point();

  Polyhedron one_ring_copy;
  typename Polyhedron::Halfedge_handle h0;
  const int nv = v->degree();
  if(nv <= MAX_DEGREE_FOR_FLIPS)
  {
    one_ring_triangulation_convex_via_flips(dual,v,one_ring_copy,h0);
  }else
  {
    // Large enough degree might be so rare that this is irrelevant.
    one_ring_triangulation_convex_via_convex_hull(dual,v,one_ring_copy,h0);
  }


  
  
  

  
  // Check that p is on the non-negative side of all faces
#ifndef NDEBUG
  for(auto f = one_ring_copy.facets_begin(); f != one_ring_copy.facets_end(); ++f)
  {
    auto h = f->halfedge();
    const auto ori = CGAL::orientation(
      h->vertex()->point(),
      h->next()->vertex()->point(),
      h->next()->next()->vertex()->point(),
      p);
    if(ori == CGAL::NEGATIVE)
    {
      throw std::runtime_error("Removed point lies on negative side of a new face, something is wrong.");
    }
  }
#endif

  auto [primal_volume, dual_volume, contains_origin] = primal_volume_subtended(one_ring_copy,p);
  //printf("primal_volume = %g, dual_volume = %g, contains_origin = %d\n",
  //  CGAL::to_double(primal_volume),
  //  CGAL::to_double(dual_volume),
  //  contains_origin);
  const Scalar cost = primal_volume;

  // Print all the triangles
  //print_faces(one_ring_copy);

  Record record;
  auto & path = record.path;
  auto & start_vertex_id = record.start_vertex_id;
  std::tie(start_vertex_id, path) = clip_ears(one_ring_copy,h0);
  //{
  //  // print path
  //  printf("start_vertex_id = %d\n",start_vertex_id);
  //  std::cout<<"path = [";
  //  for(const auto & p : path)
  //  {
  //    printf("%d ",p);
  //  }
  //  std::cout<<"];"<<std::endl;
  //}

  return {cost, record};
}

template <class Vertex_handle>
std::vector<Vertex_handle>
collect_neighbors(const Vertex_handle & v)
{
  std::vector<Vertex_handle> neighbors;
  auto h = v->halfedge();
  do
  {
    neighbors.push_back(h->opposite()->vertex());
    h = h->next()->opposite();
  } while(h != v->halfedge());
  return neighbors;
}

template <class Polyhedron, typename x0_type>
std::tuple<
    Eigen::Matrix<typename Polyhedron::Traits::FT, Eigen::Dynamic, 3, Eigen::RowMajor>,
    Eigen::VectorXi,
    Eigen::VectorXi>
dual_to_primal_mesh(
  const Polyhedron & dual,
  const x0_type & x0_exact)
{
  using Scalar = typename Polyhedron::Traits::FT;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> pV(dual.size_of_facets(), 3);
  {
    // Assign id to every face
    for(auto f = dual.facets_begin(); f != dual.facets_end(); ++f)
    {
      // get face normal and barycenter
      auto h = f->halfedge();
      const auto & A = h->vertex()->point();
      const auto & B = h->next()->vertex()->point();
      const auto & C = h->next()->next()->vertex()->point();
      const auto n = CGAL::cross_product(B - A, C - A);
      const auto bc = ((A - CGAL::ORIGIN) + (B - CGAL::ORIGIN) + (C - CGAL::ORIGIN)) / Scalar(3);
      const Scalar beta = n * bc;
      // convert x0_exact to same kernel as dual if needed
      using EK = typename CGAL::Kernel_traits<x0_type>::Kernel;
      typename Polyhedron::Traits::Point_3 x0;
      {
        CGAL::Cartesian_converter<EK, typename Polyhedron::Traits> to_dual_kernel;
        x0 = to_dual_kernel(x0_exact);
      }
      // p = x0 + n/beta
      const auto p = CGAL::ORIGIN + (x0 - CGAL::ORIGIN) + (n / beta);
      pV.row(f->id()) << p.x(), p.y(), p.z();
    }
  }
  std::vector<int> pPI;
  std::vector<int> pPC = {0};
  {
    // loop over every vertex
    for(auto v = dual.vertices_begin(); v != dual.vertices_end(); ++v)
    {
      auto h = v->halfedge();
      int np = 0;
      // loop over every face adjacent to vertex
      do
      {
        auto f = h->facet();
        pPI.push_back(f->id());
        h = h->next()->opposite();
        np++;
      } while(h != v->halfedge());
      pPC.push_back(pPC.back() + np);
    }
  }
  return {
    pV,
    Eigen::VectorXi(Eigen::VectorXi::Map(pPI.data(), pPI.size())),
    Eigen::VectorXi(Eigen::VectorXi::Map(pPC.data(), pPC.size()))};
}


int main(int argc, char *argv[])
{
  // get MAX_DEGREE_FOR_FLIPS from environment
  const char * max_degree_for_flips_env = std::getenv("MAX_DEGREE_FOR_FLIPS");
  if(max_degree_for_flips_env)
  {
    MAX_DEGREE_FOR_FLIPS = std::stoi(max_degree_for_flips_env);
  }
  printf("MAX_DEGREE_FOR_FLIPS = %d\n", MAX_DEGREE_FOR_FLIPS);

  IGL_TICTOC_LAMBDA;
  // Initial input mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(argc > 1)
  {
    igl::read_triangle_mesh(argv[1], V, F);
  }else
  {
    igl::icosahedron(V, F);
  }

  printf("|V| = %d, |F| = %d\n", V.rows(), F.rows());


  using EK = CGAL::Exact_predicates_exact_constructions_kernel;
  using EPolyhedron = CGAL::Polyhedron_3<EK, CGAL::Polyhedron_items_with_id_3>;
  using EPoint = typename EPolyhedron::Traits::Point_3;
  using IK = CGAL::Exact_predicates_inexact_constructions_kernel;
  //using IK = CGAL::Exact_predicates_exact_constructions_kernel;
  using IPolyhedron = CGAL::Polyhedron_3<IK, CGAL::Polyhedron_items_with_id_3>;
  using IPoint = typename IPolyhedron::Traits::Point_3;

  tictoc();
  auto primal_points = point_list<EK>(V);
  EPolyhedron primal;
  CGAL::convex_hull_3(primal_points.begin(),primal_points.end(),primal);
  printf("primal: %d vertices, %d facets\n",primal.size_of_vertices(),primal.size_of_facets());
  printf("initial primal convex hull: %g secs\n",tictoc());

  tictoc();
  double primal_squared_area_tol = 1e-15;
  double dual_min_edge_tol = 1e-7;
  // use chebyshev_center.h but only consider halfspaces of faces with squard areas
  // greater than primal_squared_area_tol, return as a cgal point in same kernel as
  // primal (even if we need to igl::assign to double to call chebyshev_center)
  // look at robust_dual_hull for an example
  EPoint x0_exact = chebyshev_center(primal,primal_squared_area_tol);
  // Using x0, convert ever face with squared area greater than
  // primal_squared_area_tol into a halfspace, then into a dual point about x0
  // using:
  // p = [n b]
  // d = -n / (n⋅x0 + b)
  // Conduct this all in the same kernel as primal.
  auto dual_points_exact = dual_points_list(primal,primal_squared_area_tol,x0_exact);
  auto dV_exact = point_matrix<EK>(dual_points_exact);
  printf("|dV_exact| = %d out of %d\n", dV_exact.rows(), primal.size_of_facets());

  // Convert to double as best we can.
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> dV;
  igl::copyleft::cgal::assign(dV_exact,true,dV);
  {
    Eigen::VectorXi _1,_2;
    igl::remove_duplicate_vertices(
      Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor>(dV),
      1e-14,
      dV,
      _1,_2);
  }
  printf("|dV| = %d out of %d\n", dV.rows(), dV_exact.rows());

  IPolyhedron dual;
  bool use_exact_kernel_for_dual_convex_hull = false;
  if(use_exact_kernel_for_dual_convex_hull)
  {
    // This doesn't seem to help.
    auto dual_points = point_list<EK>(dV);
    EPolyhedron dual_exact;
    CGAL::convex_hull_3(dual_points.begin(),dual_points.end(),dual_exact);
    // Copy to inexact kernel to do further processing
    {
      CGAL::Polyhedron_copy_3<EPolyhedron, IPolyhedron::HalfedgeDS> modifier(dual_exact);
      dual.delegate(modifier);
      CGAL::Cartesian_converter<EK, IK> to_inexact;
      auto v_src = dual_exact.vertices_begin();
      auto v_dst = dual.vertices_begin();
      for(; v_src != dual_exact.vertices_end(); ++v_src, ++v_dst)
      {
        v_dst->point() = to_inexact(v_src->point());
      }
    }
  }
  else
  {
    auto dual_points = point_list<IK>(dV);
    CGAL::convex_hull_3(dual_points.begin(),dual_points.end(),dual);
  }
  printf("dual convex hull: %g secs\n",tictoc());

  // wait for user to push enter, to give time to attach profiler if desired
  //primal_points.clear();
  //{
  //  std::cout<<"Press Enter to continue..."<<std::endl;
  //  std::cin.get();
  //}

  tictoc();
  // Set ids of all vertices
  {
    int vid = 0;
    for(auto v = dual.vertices_begin(); v != dual.vertices_end(); ++v)
    {
      v->id() = vid++;
    }
  }

  //igl::copyleft::cgal::join_coplanar_neighboring_facets(dual);

  //printf("confirm_all_edges_are_convex(dual)\n");
#ifndef NDEBUG
  confirm_all_edges_are_convex(dual);
#endif
  printf("dual: %d vertices, %d facets\n",dual.size_of_vertices(),dual.size_of_facets());

  using Scalar = decltype(dual)::Traits::FT;
  const int num_vertices = dual.size_of_vertices();
  struct FullRecord
  {
    // I'm skeptical whether we need to store both the latest cost _and_ a
    // "visit" counter.
    Scalar cost;
    int visit;
    typename decltype(dual)::Vertex_handle vertex;
    Record record;
  };

  //num_flips_list.clear();
  //num_flips_list.reserve(num_vertices);
  //local_size_list.clear();
  //local_size_list.reserve(num_vertices);

  std::vector<FullRecord> full_records(num_vertices);
  igl::min_heap<std::tuple<Scalar,int,int>> Q;
  for(auto v = dual.vertices_begin(); v != dual.vertices_end(); ++v)
  {
     auto [cost, record] = measure_vertex_erasure(dual, v);
     int vid = v->id();
     full_records[vid] = {cost, 0, v, record};
     Q.push({cost, vid, 0});
  }

  assert(dual.is_pure_triangle());
  assert(dual.is_closed());
  int num_vertices_removed = 0;
  const int initial_num_vertices = dual.size_of_vertices();
  const int target_num_vertices = 18;
  const int max_vertices_removed = initial_num_vertices - target_num_vertices;
  while(true)
  {
    int id;
    Scalar cost;
    int visit;
    while(true)
    {
      std::tie(cost, id, visit) = Q.top();
      Q.pop();
      if(visit == full_records[id].visit) { break; }
      //printf("skipping stale entry on %d\n", id);
    }
    if(cost == std::numeric_limits<Scalar>::infinity())
    {
      printf("Remaining entries in Q are ∞ cost\n");
      break;
    }
    //printf("removing from Q vertex %d with cost %g and visit %d\n", id, cost,visit);
    assert(cost == full_records[id].cost);
    auto v = full_records[id].vertex;
    auto record = full_records[id].record;
    const auto neighbors = collect_neighbors(v);

    //{
    //  printf("  neighbors:");
    //  for(const auto & n : neighbors)
    //  {
    //    printf(" %d", n->id());
    //  }
    //  printf("\n");
    //}

    //print_VF(dual,"b");
    erase_vertex_and_clip_ears(dual, v, record.start_vertex_id, record.path);
    num_vertices_removed++;
    //print_VF(dual,"a");
    assert(dual.is_pure_triangle());
    assert(dual.is_closed());

    for(const auto & n : neighbors)
    {
      const int nid = n->id();
      std::tie(full_records[nid].cost, full_records[nid].record) = measure_vertex_erasure(dual, n);
      full_records[nid].visit++;
      //printf("  updating neighbor %d with new cost %g and visit %d\n", nid, full_records[nid].cost, full_records[nid].visit);
      Q.push({full_records[nid].cost, nid, full_records[nid].visit});
    }

    if(num_vertices_removed >= max_vertices_removed)
    {
      printf("Removed %d vertices, stopping.\n", num_vertices_removed);
      break;
    }
  }
  printf("simplification: %g secs\n",tictoc());

  tictoc();
  // Assign id to every face
  {
    int fid = 0;
    for(auto f = dual.facets_begin(); f != dual.facets_end(); ++f)
    {
      f->id() = fid++;
    }
  }
  auto [pV,pPI,pPC] = dual_to_primal_mesh(dual,x0_exact);
  printf("dual_to_primal_mesh: %g secs\n",tictoc());
  printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  std::cout<<igl::matlab_format(pV,"pV")<<std::endl;
  std::cout<<igl::matlab_format_index(pPI.transpose().eval(),"pPI")<<std::endl;
  std::cout<<igl::matlab_format(      pPC.transpose().eval(),"pPC")<<std::endl;
  //printf("---------------------------------------------------\n");
  //std::cout<<igl::matlab_format(Eigen::VectorXi::Map(num_flips_list.data(), num_flips_list.size()).transpose().eval(),"num_flips_list")<<std::endl;
  //std::cout<<igl::matlab_format(Eigen::VectorXi::Map(local_size_list.data(), local_size_list.size()).transpose().eval(),"local_size_list")<<std::endl;
  //std::cout<<igl::matlab_format(Eigen::VectorXd::Map(flip_times_list.data(), flip_times_list.size()).transpose().eval(),"flip_times_list")<<std::endl;
  //igl::writeDMAT("num_flips_list.dmat", Eigen::VectorXi::Map(num_flips_list.data(), num_flips_list.size()));
  //igl::writeDMAT("local_size_list.dmat", Eigen::VectorXi::Map(local_size_list.data(), local_size_list.size()));
  //igl::writeDMAT("flip_times_list.dmat", Eigen::VectorXd::Map(flip_times_list.data(), flip_times_list.size()));


}
