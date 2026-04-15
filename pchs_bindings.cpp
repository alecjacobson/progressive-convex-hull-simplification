#include "convex_hull_simplification.h"

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(pchs, m)
{
  m.doc() = "Progressive Convex Hull Simplification";

  // --- CostFunction enum ---
  nb::enum_<CostFunction>(m, "CostFunction")
    .value("volume", CostFunction::PRIMAL_VOLUME,     "Primal volume added per removal")
    .value("area",   CostFunction::PRIMAL_AREA_CHANGE, "Change in primal surface area per removal")
    .export_values();

  // --- Stats struct ---
  nb::class_<ConvexHullSimplification::Stats>(m, "Stats")
    .def_ro("t_primal_hull",   &ConvexHullSimplification::Stats::t_primal_hull,
            "Time to compute the initial primal convex hull (s)")
    .def_ro("t_dual_hull",     &ConvexHullSimplification::Stats::t_dual_hull,
            "Time to compute the dual convex hull (s)")
    .def_ro("t_queue_init",    &ConvexHullSimplification::Stats::t_queue_init,
            "Time to initialize the priority queue (s)")
    .def_ro("t_last_simplify", &ConvexHullSimplification::Stats::t_last_simplify,
            "Time spent in the most recent simplify_to() call (s)");

  // --- ConvexHullSimplification class ---
  nb::class_<ConvexHullSimplification>(m, "ConvexHullSimplification",
    "Progressively simplifies a convex hull by greedily removing dual vertices.\n\n"
    "Example::\n\n"
    "  chs = pchs.ConvexHullSimplification(V, F)\n"
    "  chs.simplify_to(18)\n"
    "  pV, pPI, pPC = chs.get_primal_mesh()\n")
    .def(nb::init<
           const Eigen::MatrixXd &,
           const Eigen::MatrixXi &,
           int,
           CostFunction>(),
         "V"_a, "F"_a,
         "max_degree_for_flips"_a = 100,
         "cost_function"_a = CostFunction::PRIMAL_VOLUME,
         "Build the primal hull, compute the Chebyshev center, construct the dual,\n"
         "and initialize the priority queue.")
    .def("step",
         &ConvexHullSimplification::step,
         "Pop and apply the cheapest vertex erasure. Returns False when exhausted.")
    .def("simplify_to",
         &ConvexHullSimplification::simplify_to,
         "target"_a,
         "Call step() until num_dual_vertices() == target or step() returns False.")
    .def("num_dual_vertices",
         &ConvexHullSimplification::num_dual_vertices,
         "Number of dual vertices (halfspaces) currently remaining.")
    .def("dual_vertex_ids",
         &ConvexHullSimplification::dual_vertex_ids,
         "Stable vertex ids in the same order as the polygons from get_primal_mesh().")
    .def("dual_greedy_coloring",
         &ConvexHullSimplification::dual_greedy_coloring,
         "k"_a = 9,
         "Greedy graph coloring of the dual 1-skeleton with at most k colors.\n"
         "Returns a VectorXi of size = initial vertex count; removed vertices get -1.")
    .def("get_primal_mesh",
         &ConvexHullSimplification::get_primal_mesh,
         "Reconstruct the primal polygon mesh as (pV, pPI, pPC).\n\n"
         "pV  — #dual_facets x 3 primal vertex positions\n"
         "pPI — flat polygon vertex index list\n"
         "pPC — cumulative polygon sizes (len = num_dual_vertices + 1)")
    .def("get_dual_mesh",
         &ConvexHullSimplification::get_dual_mesh,
         "Extract the current dual as a triangle mesh (V, F) with remapped indices.")
    .def("popped_dual_vertex_costs",
         &ConvexHullSimplification::popped_dual_vertex_costs,
         "Cost at which each vertex was popped; NaN for vertices still alive.\n"
         "Indexed by original vertex id (size = initial vertex count).")
    .def("popped_dual_vertex_ids",
         &ConvexHullSimplification::popped_dual_vertex_ids,
         "Vertex ids in the order they were removed from the dual.")
    .def("stats",
         &ConvexHullSimplification::stats,
         "Timing breakdown for each construction/simplification phase.");

  // --- Convenience free function ---
  m.def("simplify_convex_hull",
    &simplify_convex_hull,
    "V"_a, "F"_a, "target"_a,
    "max_degree_for_flips"_a = 100,
    "cost_function"_a = CostFunction::PRIMAL_VOLUME,
    "Construct, simplify to target dual vertices, and return (pV, pPI, pPC).");
}
