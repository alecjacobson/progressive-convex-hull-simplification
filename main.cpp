#include "convex_hull_simplification.h"

#include <igl/read_triangle_mesh.h>
#include <igl/icosahedron.h>
#include <igl/matlab_format.h>
#include <igl/get_seconds.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/assign.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>
#include <string>
#include <cstdio>
#include <cstdlib>

// --- Debug helpers (all call sites are commented out) ---

template <class Polyhedron>
void print_in_matlab_format(const Polyhedron & poly, const std::string prefix = "")
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
    printf("v->id() = %zu\n",v->id());
    if(int(v->id()) > max_id) { max_id = int(v->id()); }
  }
  printf("-----------------\n");
  printf("max_id = %d\n",max_id);
  Eigen::Matrix<typename Polyhedron::Traits::FT,Eigen::Dynamic,3,Eigen::RowMajor> V(max_id+1,3);
  V.setConstant(std::numeric_limits<typename Polyhedron::Traits::FT>::quiet_NaN());
  for(auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    V.row(v->id()) << v->point().x(), v->point().y(), v->point().z();
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
    assert(k == (int)poly.size_of_facets());
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
    printf("  (%zu %zu %zu)\n",
      h->vertex()->id(),
      h->next()->vertex()->id(),
      h->next()->next()->vertex()->id());
  }
}

// --- Polygon mesh triangulation ---

// Fan-triangulate a polygon mesh given in (pV, pPI, pPC) format.
// Each polygon p spans pPI[pPC[p] .. pPC[p+1]-1].
std::pair<Eigen::MatrixXd, Eigen::MatrixXi>
triangulate_polygon_mesh(
  const Eigen::MatrixXd & pV,
  const Eigen::VectorXi & pPI,
  const Eigen::VectorXi & pPC)
{
  int nf = 0;
  for(int p = 0; p < pPC.size()-1; p++)
    nf += pPC(p+1) - pPC(p) - 2;

  Eigen::MatrixXi F(nf, 3);
  int fi = 0;
  for(int p = 0; p < pPC.size()-1; p++)
  {
    const int first = pPI(pPC(p));
    for(int j = pPC(p)+1; j < pPC(p+1)-1; j++)
      F.row(fi++) << first, pPI(j), pPI(j+1);
  }
  return {pV, F};
}

// --- Entry point ---

int main(int argc, char *argv[])
{
  int max_degree_for_flips = 100;
  if(const char * env = std::getenv("MAX_DEGREE_FOR_FLIPS"))
    max_degree_for_flips = std::stoi(env);

  // Parse flags
  bool print_stats = false;
  bool interactive = false;
  int target = 18;
  const char * mesh_file = nullptr;
  for(int i = 1; i < argc; i++)
  {
    const std::string arg(argv[i]);
    if(arg == "--stats")       { print_stats = true; }
    else if(arg == "--interactive") { interactive = true; }
    else if(arg == "--target" && i+1 < argc) { target = std::stoi(argv[++i]); }
    else                       { mesh_file = argv[i]; }
  }

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(mesh_file)
    igl::read_triangle_mesh(mesh_file, V, F);
  else
    igl::icosahedron(V, F);

  printf("|V| = %d, |F| = %d\n", V.rows(), F.rows());

  // --- Interactive mode ---
  if(interactive)
  {
    using Viewer = igl::opengl::glfw::Viewer;

    auto chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips);

    Viewer viewer;

    // Mesh 0: input mesh
    viewer.data(0).set_mesh(V, F);
    viewer.data(0).set_colors((Eigen::MatrixXd(1,3) << 0.7, 0.75, 0.8).finished());

    // Mesh 1: simplified hull
    const int hull_id = static_cast<int>(viewer.append_mesh());

    // Helper: push current hull state into the viewer
    const auto update_hull = [&]()
    {
      auto [pV, pPI, pPC] = chs->get_primal_mesh();
      auto [hV, hF] = triangulate_polygon_mesh(pV, pPI, pPC);
      viewer.data(hull_id).clear();
      viewer.data(hull_id).set_mesh(hV, hF);
      // RGBA: warm orange at 75% opacity
      viewer.data(hull_id).set_colors(
        (Eigen::MatrixXd(1,4) << 1.0, 0.6, 0.2, 0.75).finished());
      printf("dual vertices: %d\n", chs->num_dual_vertices());
    };
    update_hull();

    viewer.callback_key_pressed =
      [&](Viewer &, unsigned char key, int) -> bool
    {
      switch(key)
      {
        case ' ':
          chs->step();
          update_hull();
          return true;
        case 'a': case 'A':
          chs->simplify_to(target);
          update_hull();
          return true;
        case 'r': case 'R':
          chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips);
          update_hull();
          return true;
      }
      return false;
    };

    printf("Keys: [space] step  [a] simplify_to(%d)  [r] reset\n", target);
    viewer.launch();
    return 0;
  }

  // --- Non-interactive mode ---
  ConvexHullSimplification chs(V, F, max_degree_for_flips);
  chs.simplify_to(target);

  if(print_stats)
  {
    const auto s = chs.stats();
    printf("primal hull:  %8.4f s\n", s.t_primal_hull);
    printf("dual hull:    %8.4f s\n", s.t_dual_hull);
    printf("queue init:   %8.4f s\n", s.t_queue_init);
    printf("simplify_to:  %8.4f s\n", s.t_last_simplify);
    printf("total:        %8.4f s\n",
      s.t_primal_hull + s.t_dual_hull + s.t_queue_init + s.t_last_simplify);
  }

  auto [pV, pPI, pPC] = chs.get_primal_mesh();
  printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  std::cout<<igl::matlab_format(pV,"pV")<<std::endl;
  std::cout<<igl::matlab_format_index(pPI.transpose().eval(),"pPI")<<std::endl;
  std::cout<<igl::matlab_format(      pPC.transpose().eval(),"pPC")<<std::endl;
}
