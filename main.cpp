#include "convex_hull_simplification.h"

#include <igl/read_triangle_mesh.h>
#include <igl/icosahedron.h>
#include <igl/matlab_format.h>
#include <igl/get_seconds.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/assign.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <imgui.h>

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

  // --- Interactive mode (polyscope) ---
  if(interactive)
  {
    auto chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips);

    polyscope::init();
    polyscope::options::programName = "Progressive Convex Hull Simplification";
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::giveFocusOnShow = true;

    // Input mesh: gray, no edges, smooth shading
    {
      auto* m = polyscope::registerSurfaceMesh("input mesh", V, F);
      m->setEdgeWidth(0.0);
      m->setSurfaceColor(glm::vec3(0.7f, 0.75f, 0.8f));
    }

    float hull_alpha = 0.65f;  // opacity: 1=opaque, 0=transparent

    // Helper: rebuild the simplified-hull mesh and polygon-edge curve network
    // in polyscope from the current chs state.
    // Faces are pseudocolored by stable dual-vertex id so colors are consistent
    // across simplification steps. Edges follow polygon boundaries (pPI/pPC),
    // not the fan-triangulation diagonals.
    const auto update_hull = [&]()
    {
      auto [pV, pPI, pPC] = chs->get_primal_mesh();
      auto [hV, hF] = triangulate_polygon_mesh(pV, pPI, pPC);

      // Greedy coloring → ColorBrewer Set1 (9 qualitative colors).
      // Labels wrap with % 9; in practice a convex polyhedron needs ≤ 4.
      static const double set1[9][3] = {
        {0.894, 0.102, 0.110},  // red      #E41A1C
        {0.216, 0.494, 0.722},  // blue     #377EB8
        {0.302, 0.686, 0.290},  // green    #4DAF4A
        {0.596, 0.306, 0.639},  // purple   #984EA3
        {1.000, 0.498, 0.000},  // orange   #FF7F00
        {1.000, 1.000, 0.200},  // yellow   #FFFF33
        {0.651, 0.337, 0.157},  // brown    #A65628
        {0.969, 0.506, 0.749},  // pink     #F781BF
        {0.600, 0.600, 0.600},  // gray     #999999
      };
      const auto vids     = chs->dual_vertex_ids();
      const auto coloring = chs->dual_greedy_coloring();
      int nf = 0;
      for(int p = 0; p < pPC.size()-1; p++)
        nf += pPC(p+1) - pPC(p) - 2;
      Eigen::MatrixXd tri_colors(nf, 3);
      {
        int fi = 0;
        for(int p = 0; p < pPC.size()-1; p++)
        {
          const double* c = set1[coloring(vids(p)) % 9];
          for(int j = pPC(p)+1; j < pPC(p+1)-1; j++)
            tri_colors.row(fi++) << c[0], c[1], c[2];
        }
      }

      auto* m = polyscope::registerSurfaceMesh("simplified hull", hV, hF);
      m->setTransparency(hull_alpha);
      m->setSmoothShade(false);   // per-face normals
      m->setEdgeWidth(0.0);       // edges come from the curve network below

      auto* q = m->addFaceColorQuantity("color", tri_colors);
      q->setEnabled(true);

      // Polygon boundary edges: one edge per consecutive pair in each polygon
      // (including the closing edge), built directly from pPI/pPC.
      const int ne = (int)pPI.size();  // sum of polygon degrees = total edges
      Eigen::MatrixXi edge_mat(ne, 2);
      {
        int ei = 0;
        for(int p = 0; p < pPC.size()-1; p++)
        {
          const int start = pPC(p), np = pPC(p+1) - pPC(p);
          for(int i = 0; i < np; i++)
            edge_mat.row(ei++) << pPI(start + i), pPI(start + (i+1) % np);
        }
      }
      auto* cn = polyscope::registerCurveNetwork("hull edges", pV, edge_mat);
      cn->setRadius(0.001, true);
      cn->setColor(glm::vec3(0.15f, 0.15f, 0.15f));

      printf("dual vertices: %d\n", chs->num_dual_vertices());
    };
    update_hull();

    bool animate = false;

    polyscope::state::userCallback = [&]()
    {
      // Animate: one simplification step per rendered frame.
      if(animate)
      {
        if(chs->num_dual_vertices() > target)
        {
          chs->step();
          update_hull();
        }
        if(chs->num_dual_vertices() <= target)
          animate = false;
      }

      ImGui::Text("Dual vertices: %d", chs->num_dual_vertices());
      ImGui::SetNextItemWidth(120.0f);
      ImGui::InputInt("Target", &target);
      ImGui::Checkbox("Animate to target", &animate);
      ImGui::SetNextItemWidth(120.0f);
      if(ImGui::SliderFloat("Hull alpha", &hull_alpha, 0.0f, 1.0f))
        if(auto* m = polyscope::getSurfaceMesh("simplified hull"))
          m->setTransparency(hull_alpha);

      if(ImGui::Button("Step"))
      {
        chs->step();
        update_hull();
      }
      ImGui::SameLine();
      if(ImGui::Button("Go to target"))
      {
        animate = false;
        chs->simplify_to(target);
        update_hull();
      }
      ImGui::SameLine();
      if(ImGui::Button("Reset"))
      {
        animate = false;
        chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips);
        update_hull();
      }
    };

    polyscope::show();
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
