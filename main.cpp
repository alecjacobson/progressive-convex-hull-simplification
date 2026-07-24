#include "convex_hull_simplification.h"

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/writeDMAT.h>
#include <igl/icosahedron.h>
#include <igl/matlab_format.h>
#include <igl/get_seconds.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/polygons_to_triangles.h>
#include <igl/doublearea.h>
#include <igl/centroid.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/assign.h>

#include "write_polygon_ply.h"

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
  const char * primal_output  = nullptr;
  const char * dual_output    = nullptr;
  const char * primal_initial = nullptr;
  const char * dual_initial   = nullptr;
  bool write_costs = false;
  const char * costs_output   = nullptr;  // nullptr = use default name
  bool write_popped_ids = false;
  const char * popped_ids_output = nullptr;  // nullptr = use default name
  CostFunction cost_function = CostFunction::PRIMAL_VOLUME;
  for(int i = 1; i < argc; i++)
  {
    const std::string arg(argv[i]);
    if(arg == "--help" || arg == "-h")
    {
      printf(
        "Usage: pchs [options] [mesh.ply]\n"
        "\n"
        "Defaults to an icosahedron if no mesh is given.\n"
        "\n"
        "Options:\n"
        "  --target N              Simplify to N dual vertices / halfspaces (default: 18);\n"
        "                          if N < 0, subtract from the initial count (e.g. -1 removes one face)\n"
        "  --cost-function F       Cost metric: 'volume' or 'area' (default: volume)\n"
        "  --primal-output file    Write simplified hull as a polygon PLY\n"
        "  --dual-output file      Write simplified dual as a triangle mesh\n"
        "  --primal-initial file   Write initial (unsimplified) primal hull\n"
        "  --dual-initial file     Write initial (unsimplified) dual hull\n"
        "  --costs [file.dmat]     Write per-vertex removal costs (default: <stem>-costs.dmat)\n"
        "  --popped-ids [file.dmat]\n"
        "                          Write removal order (default: <stem>-popped-ids.dmat)\n"
        "  --stats                 Print per-phase timing after simplification\n"
#ifdef PCHS_INTERACTIVE
        "  --interactive           Open a polyscope viewer\n"
#endif
        "  --help, -h              Show this message\n"
      );
      return 0;
    }
    else if(arg == "--stats")       { print_stats = true; }
    else if(arg == "--interactive")
    {
#ifdef PCHS_INTERACTIVE
      interactive = true;
#else
      fprintf(stderr, "error: --interactive requires a build with PCHS_INTERACTIVE=ON\n");
      return 1;
#endif
    }
    else if(arg == "--target"          && i+1 < argc) { target          = std::stoi(argv[++i]); }
    else if(arg == "--primal-output"   && i+1 < argc) { primal_output   = argv[++i]; }
    else if(arg == "--dual-output"     && i+1 < argc) { dual_output     = argv[++i]; }
    else if(arg == "--primal-initial"  && i+1 < argc) { primal_initial  = argv[++i]; }
    else if(arg == "--dual-initial"    && i+1 < argc) { dual_initial    = argv[++i]; }
    else if(arg == "--costs")
    {
      write_costs = true;
      if(i+1 < argc && argv[i+1][0] != '-') costs_output = argv[++i];
    }
    else if(arg == "--popped-ids")
    {
      write_popped_ids = true;
      if(i+1 < argc && argv[i+1][0] != '-') popped_ids_output = argv[++i];
    }
    else if(arg == "--cost-function" && i+1 < argc)
    {
      const std::string cf(argv[++i]);
      if(cf == "volume")     cost_function = CostFunction::PRIMAL_VOLUME;
      else if(cf == "area")  cost_function = CostFunction::PRIMAL_AREA;
      else if(cf == "mean-width")  cost_function = CostFunction::PRIMAL_MEAN_WIDTH;
      else { fprintf(stderr, "error: unknown cost function '%s' (use 'volume', 'area' or 'mean-width')\n", cf.c_str()); return 1; }
    }
    else                       { mesh_file = argv[i]; }
  }

  // Derive default output filenames from the input stem, or use fallbacks.
  auto file_stem = [](const char* path) -> std::string
  {
    std::string s(path);
    auto sep = s.find_last_of("/\\");
    if(sep != std::string::npos) s = s.substr(sep + 1);
    auto dot = s.rfind('.');
    if(dot != std::string::npos) s = s.substr(0, dot);
    return s;
  };
  const std::string stem = mesh_file ? file_stem(mesh_file) : "";
  const std::string primal_path = primal_output ? primal_output
                                : (stem.empty() ? "primal-output.ply" : (stem + "-primal.ply"));
  const std::string dual_path   = dual_output   ? dual_output
                                : (stem.empty() ? "dual-output.ply"   : (stem + "-dual.ply"));
  const std::string costs_path      = costs_output      ? costs_output
                                    : (stem.empty() ? "costs.dmat"      : (stem + "-costs.dmat"));
  const std::string popped_ids_path = popped_ids_output ? popped_ids_output
                                    : (stem.empty() ? "popped-ids.dmat" : (stem + "-popped-ids.dmat"));

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(mesh_file)
  {
    igl::read_triangle_mesh(mesh_file, V, F);
    if(V.rows() == 3*F.rows())
    {
      const int nv0 = V.rows();
      Eigen::VectorXi _1, _2;
      igl::remove_duplicate_vertices(
          Eigen::MatrixXd(V),Eigen::MatrixXi(F),0,
          V,_1,_2,F);
      printf("Merged %d duplicate vertices from STL file\n", nv0 - V.rows());
    }
  }
  else
  {
    igl::icosahedron(V, F);
  }

  printf("|V| = %d, |F| = %d\n", V.rows(), F.rows());

#ifdef PCHS_INTERACTIVE
  // --- Interactive mode (polyscope) ---
  if(interactive)
  {
    auto chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips, cost_function);
    if(target < 0) target = chs->num_dual_vertices() + target;

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
    // cost_function_idx mirrors cost_function: 0=volume, 1=area, 2=mean-width
    int cost_function_idx = (cost_function == CostFunction::PRIMAL_AREA) ? 1
                          : (cost_function == CostFunction::PRIMAL_MEAN_WIDTH)  ? 2 : 0;

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
          static const char* cf_modes[] = {"Volume", "Area", "Mean Width"};
          ImGui::SetNextItemWidth(120.0f);
          if(ImGui::Combo("Cost function", &cost_function_idx, cf_modes, 3))
          {
            animate = false;
            const int current_n = chs->num_dual_vertices();
            cost_function = (cost_function_idx == 1) ? CostFunction::PRIMAL_AREA
                          : (cost_function_idx == 2) ? CostFunction::PRIMAL_MEAN_WIDTH
                          :                            CostFunction::PRIMAL_VOLUME;
            chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips, cost_function);
            chs->simplify_to(current_n);
            update_hull();
          }
        }
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
          chs = std::make_unique<ConvexHullSimplification>(V, F, max_degree_for_flips, cost_function);
          update_hull();
        }
      }
    };

    polyscope::show();
    ao_thread.join();
    return 0;
  }
#endif // PCHS_INTERACTIVE

  // Helper: write a file and print a one-line status, or print an error.
  auto write_poly = [](const std::string & path,
                       const Eigen::MatrixXd & V,
                       const Eigen::VectorXi & PI,
                       const Eigen::VectorXi & PC)
  {
    if(!write_polygon_ply(path, V, PI, PC))
      fprintf(stderr, "error: failed to write %s\n", path.c_str());
    else
      printf("wrote: %s\n", path.c_str());
  };
  auto write_tri = [](const std::string & path,
                      const Eigen::MatrixXd & V,
                      const Eigen::MatrixXi & F)
  {
    if(!igl::write_triangle_mesh(path, V, F))
      fprintf(stderr, "error: failed to write %s\n", path.c_str());
    else
      printf("wrote: %s\n", path.c_str());
  };

  //// Wait for user to press Enter before starting simplification, to give them a
  //// chance to read the initial stats and outputs.
  //{
  //  printf("Press Enter to start simplification...");
  //  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //}

  // --- Non-interactive mode ---
  ConvexHullSimplification chs(V, F, max_degree_for_flips, cost_function);
  if(target < 0) target = chs.num_dual_vertices() + target;

  // Initial outputs (before any simplification).
  double area_initial = 0.0;
  double volume_initial = 0.0;
  double mean_width_initial = 0.0;
  const auto primal_area_and_volume = [](auto & chs) -> std::pair<double, double>
  {
    auto [pV, pPI, pPC] = chs.get_primal_mesh();
    Eigen::MatrixXi pF;
    Eigen::VectorXi _;
    igl::polygons_to_triangles(pPI, pPC, pF, _);
    Eigen::VectorXd pA;
    igl::doublearea(pV, pF, pA);
    double area = 0.5 * pA.sum();
    double volume;
    {
      Eigen::RowVector3d _;
      igl::centroid(pV,pF,_,volume);
    }
    return {area, volume};
  };
  if(print_stats)
  {
    std::tie(area_initial, volume_initial) = primal_area_and_volume(chs);
    mean_width_initial = chs.mean_width();
  }

  if(primal_initial)
  {
    auto [ipV, ipPI, ipPC] = chs.get_primal_mesh();
    write_poly(primal_initial, ipV, ipPI, ipPC);
  }
  
  if(dual_initial)
  {
    auto [idV, idF] = chs.get_dual_mesh();
    write_tri(dual_initial, idV, idF);
  }

  chs.simplify_to(target);

  if(print_stats)
  {
    auto [area_final, volume_final] = primal_area_and_volume(chs);
    const double mean_width_final = chs.mean_width();
    const Eigen::VectorXd costs = chs.popped_dual_vertex_costs();
    const Eigen::VectorXi ids = chs.popped_dual_vertex_ids();
    //printf("last id: %d\n",ids(ids.size()-1));

    const auto s = chs.stats();
    printf("primal hull:   %8.4g s\n", s.t_primal_hull);
    printf("dual hull:     %8.4g s\n", s.t_dual_hull);
    printf("queue init:    %8.4g s\n", s.t_queue_init);
    printf("simplify_to:   %8.4g s\n", s.t_last_simplify);
    printf("total:         %8.4g s\n",
      s.t_primal_hull + s.t_dual_hull + s.t_queue_init + s.t_last_simplify);
    printf("initial area:  %8.4g\n", area_initial);
    printf("final area:    %8.4g\n", area_final);
    if(cost_function == CostFunction::PRIMAL_AREA)
    {
      printf("∫ final area:  %8.4g\n", area_initial + costs(ids).sum());
    }
    printf("f/i area:      %8.4g\n", area_final / area_initial);
    printf("initial volume:%8.4g\n", volume_initial);
    printf("final volume:  %8.4g\n", volume_final);
    if(cost_function == CostFunction::PRIMAL_VOLUME)
    {
      printf("∫ vol change:  %8.4g\n", volume_initial + costs(ids).sum());
    }
    printf("f/i volume:    %8.4g\n", volume_final / volume_initial);
    printf("initial mw:    %8.4g\n", mean_width_initial);
    printf("final mw:      %8.4g\n", mean_width_final);
    if(cost_function == CostFunction::PRIMAL_MEAN_WIDTH)
    {
      printf("∫ mw change:   %8.4g\n", mean_width_initial + costs(ids).sum());
    }
    printf("f/i mw:        %8.4g\n", mean_width_final / mean_width_initial);
  }

  if(write_costs)
  {
    const Eigen::VectorXd costs = chs.popped_dual_vertex_costs();
    if(!igl::writeDMAT(costs_path, costs))
      fprintf(stderr, "error: failed to write %s\n", costs_path.c_str());
    else
      printf("wrote: %s\n", costs_path.c_str());
  }

  if(write_popped_ids)
  {
    const Eigen::VectorXi ids = chs.popped_dual_vertex_ids();
    if(!igl::writeDMAT(popped_ids_path, ids))
      fprintf(stderr, "error: failed to write %s\n", popped_ids_path.c_str());
    else
      printf("wrote: %s\n", popped_ids_path.c_str());
  }

  {
    auto [pV, pPI, pPC] = chs.get_primal_mesh();
    write_poly(primal_path, pV, pPI, pPC);
  }
  {
    auto [dV, dF] = chs.get_dual_mesh();
    write_tri(dual_path, dV, dF);
  }
}
