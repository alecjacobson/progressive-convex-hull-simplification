#include "convex_hull_simplification.h"

#include <igl/read_triangle_mesh.h>
#include <igl/icosahedron.h>
#include <igl/matlab_format.h>
#include <igl/get_seconds.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/assign.h>

#ifdef PCHS_INTERACTIVE
#include <igl/polygons_to_triangles.h>
#include <igl/per_vertex_normals.h>
#include <igl/embree/ambient_occlusion.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <imgui.h>
#endif // PCHS_INTERACTIVE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <memory>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#ifdef PCHS_INTERACTIVE
#include <atomic>
#include <mutex>
#include <thread>
#endif // PCHS_INTERACTIVE

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
    else if(arg == "--interactive")
    {
#ifdef PCHS_INTERACTIVE
      interactive = true;
#else
      fprintf(stderr, "error: --interactive requires a build with PCHS_INTERACTIVE=ON\n");
      return 1;
#endif
    }
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

#ifdef PCHS_INTERACTIVE
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
      m->setEnabled(true);
    }

    // AO: computed in a worker thread; results picked up in userCallback.
    std::atomic<bool> ao_ready{false};
    Eigen::VectorXd   ao_values;
    std::mutex        ao_mutex;
    std::thread ao_thread([&]()
    {
      Eigen::MatrixXd N;
      igl::per_vertex_normals(V, F, N);
      Eigen::VectorXd S;
      igl::embree::ambient_occlusion(V, F, V, N, 512, S);
      {
        std::lock_guard<std::mutex> lk(ao_mutex);
        ao_values = (1.0 - S.array()).matrix();  // 1 = unoccluded (bright)
      }
      ao_ready.store(true, std::memory_order_release);
    });

    float hull_alpha = 0.85f;  // opacity: 1=opaque, 0=transparent
    int   color_mode = 1;      // 0 = graph coloring, 1 = normal pseudocolor

    // Helper: rebuild the simplified-hull mesh and polygon-edge curve network
    // in polyscope from the current chs state.
    const auto update_hull = [&]()
    {
      auto [pV, pPI, pPC] = chs->get_primal_mesh();
      Eigen::MatrixXf pF;
      Eigen::VectorXi pJ;
      igl::polygons_to_triangles(pPI, pPC, pF, pJ);

      Eigen::MatrixXd tri_colors(pF.rows(),3);

      if(color_mode == 0)
      {
        // ColorBrewer Set1 graph coloring (9 qualitative colors).
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
        const auto coloring = chs->dual_greedy_coloring(9);
        int fi = 0;
        for(int p = 0; p < pPC.size()-1; p++)
        {
          const double* c = set1[coloring(vids(p)) % 9];
          for(int j = pPC(p)+1; j < pPC(p+1)-1; j++)
            tri_colors.row(fi++) << c[0], c[1], c[2];
        }
      }
      else
      {
        Eigen::MatrixXd pN = Eigen::MatrixXd::Zero(pPC.size()-1, 3);
        for(int f = 0; f < pF.rows(); f++)
        {
          // area-weighted normal for this face
          const auto nf = (pV.row(pF(f,1)) - pV.row(pF(f,0))).head<3>().cross(
                          (pV.row(pF(f,2)) - pV.row(pF(f,0))).head<3>()).eval();
          pN.row(pJ(f)) += nf;
        }
        for(int i = 0; i < pN.rows(); i++)
        {
          pN.row(i) = (0.5*pN.row(i).normalized().array() + 0.5).eval();
        }
        for(int f = 0; f < pF.rows(); f++)
        {
          tri_colors.row(f) = pN.row(pJ(f));
        }
      }

      auto* m = polyscope::registerSurfaceMesh("simplified hull", pV, pF);
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
      // On the first frame, collapse polyscope's built-in menus.
      // Uses ImGui state storage; takes effect from frame 2.
      static bool first_frame = true;
      if(first_frame)
      {
        first_frame = false;
        ImGui::GetStateStorage()->SetInt(ImGui::GetID("Polyscope"),  0);
        ImGui::GetStateStorage()->SetInt(ImGui::GetID("Structures"), 0);
      }

      // Pick up AO result as soon as the worker thread finishes.
      if(ao_ready.load(std::memory_order_acquire))
      {
        ao_ready.store(false, std::memory_order_relaxed);
        Eigen::VectorXd S;
        {
          std::lock_guard<std::mutex> lk(ao_mutex);
          S = ao_values;
        }
        if(auto* m = polyscope::getSurfaceMesh("input mesh"))
        {
          const glm::vec3 base(0.7f, 0.75f, 0.8f);
          Eigen::MatrixXd C(S.size(), 3);
          for(int i = 0; i < S.size(); i++)
            C.row(i) << base.r * S(i), base.g * S(i), base.b * S(i);
          auto* q = m->addVertexColorQuantity("AO", C);
          q->setEnabled(true);
        }
      }

      // Animation runs unconditionally so it continues even when the panel
      // is collapsed. ceil(remaining^0.75) steps per frame.
      if(animate)
      {
        const int remaining = chs->num_dual_vertices() - target;
        if(remaining > 0)
        {
          const int n = std::max(1, (int)std::ceil(std::pow((double)remaining, 0.75)));
          for(int s = 0; s < n && chs->num_dual_vertices() > target; ++s)
            chs->step();
          update_hull();
        }
        if(chs->num_dual_vertices() <= target)
          animate = false;
      }

      if(ImGui::CollapsingHeader("Hull Simplification", ImGuiTreeNodeFlags_DefaultOpen))
      {
        ImGui::Text("Dual vertices: %d", chs->num_dual_vertices());
        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputInt("Target", &target);
        ImGui::Checkbox("Animate to target", &animate);
        ImGui::SetNextItemWidth(120.0f);
        if(ImGui::SliderFloat("Hull alpha", &hull_alpha, 0.0f, 1.0f))
          if(auto* m = polyscope::getSurfaceMesh("simplified hull"))
            m->setTransparency(hull_alpha);
        {
          static const char* modes[] = {"Graph coloring", "Normal pseudocolor"};
          ImGui::SetNextItemWidth(180.0f);
          if(ImGui::Combo("Coloring", &color_mode, modes, 2))
            update_hull();
        }

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
      }
    };

    polyscope::show();
    ao_thread.join();
    return 0;
  }
#endif // PCHS_INTERACTIVE

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
