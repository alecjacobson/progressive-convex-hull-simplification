#pragma once
#include "vertex_erasure.h"
#include "mesh_backend.h"

#include <igl/min_heap.h>

#include <Eigen/Core>
#include <tuple>
#include <vector>

// Progressively simplifies a convex hull by greedily removing the dual vertex
// (= halfspace) that adds the least primal volume at each step.
//
// Usage:
//   ConvexHullSimplification chs(V, F);
//   while(chs.step()) { auto [pV,pPI,pPC] = chs.get_primal_mesh(); ... }
//   chs.simplify_to(18);
//   auto [pV,pPI,pPC] = chs.get_primal_mesh();
class ConvexHullSimplification
{
public:
  // Kernel/mesh types come from the selected backend (see mesh_backend.h).
  using EK = mesh::Kernel;
  using IK = mesh::Kernel;
  using IPolyhedron = mesh::Polyhedron;
  using Scalar = mesh::Kernel::FT;

  // Build the primal convex hull, compute the Chebyshev center, construct the
  // dual, and initialize the priority queue with costs for every dual vertex.
  ConvexHullSimplification(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    int max_degree_for_flips = 100,
    CostFunction cost_function = CostFunction::PRIMAL_VOLUME);

  // Drain stale queue entries, then pop and apply the cheapest non-stale
  // vertex erasure. Updates neighbor costs. Returns false if the queue is
  // exhausted or every remaining vertex has infinite cost.
  bool step();

  // Repeatedly call step() until num_dual_vertices() == target or step()
  // returns false. No-op if already at or below target.
  void simplify_to(int target_num_dual_vertices);

  // Number of dual vertices currently remaining (= number of halfspaces).
  int num_dual_vertices() const;

  // Original vertex ids of the remaining dual vertices, in the same order as
  // the polygons produced by get_primal_mesh() (i.e. polygon i has id result[i]).
  // Ids are stable across simplification steps: a surviving vertex always keeps
  // the id it was assigned at construction.
  Eigen::VectorXi dual_vertex_ids() const;

  // Greedy graph coloring of the current dual: assigns an integer label in
  // [0, k) to each dual vertex such that no two adjacent vertices share a
  // label. Returns a VectorXi of size = initial vertex count, indexed by
  // vertex id; removed vertices are assigned -1. k defaults to 9.
  Eigen::VectorXi dual_greedy_coloring(int k = 9) const;

  // Reconstruct the primal polygon mesh from the current dual state.
  // pV  — #dual_facets x 3 primal vertex positions (one per dual face)
  // pPI — flat polygon vertex index list
  // pPC — cumulative counts (size = num_dual_vertices() + 1)
  std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi> get_primal_mesh();

  // Extract the current dual mesh as a triangle mesh (V, F).
  // Vertex ids are remapped to a contiguous 0-based range.
  std::pair<Eigen::MatrixXd, Eigen::MatrixXi> get_dual_mesh() const;

  Eigen::VectorXd popped_dual_vertex_costs() const;
  Eigen::VectorXi popped_dual_vertex_ids() const;

  struct Stats
  {
    double t_primal_hull;    // convex_hull_3 on input points
    double t_dual_hull;      // Chebyshev center + dual points + dedup + convex_hull_3
    double t_queue_init;     // initial measure_vertex_erasure for all dual vertices
    double t_last_simplify;  // most recent simplify_to() call (0 if never called)
  };
  Stats stats() const;

  // Mean width of the current primal hull: (1/4π) Σ ℓ_e θ_e over primal edges,
  // where θ_e is the exterior dihedral angle. Recomputed from scratch each call.
  double mean_width() const;

private:
  struct FullRecord
  {
    Scalar cost;
    int visit;
    IPolyhedron::Vertex_handle vertex;
    Record record;
  };

  EK::Point_3 x0_exact_;
  IPolyhedron dual_;
  std::vector<FullRecord> full_records_;
  std::vector<int> pop_ids_;
  igl::min_heap<std::tuple<Scalar,int,int>> Q_;
  int max_degree_for_flips_;
  CostFunction cost_function_;
  Stats stats_;
};

// Convenience wrapper: construct, simplify to target, return polygon mesh.
std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi>
simplify_convex_hull(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int target_num_dual_vertices,
  int max_degree_for_flips = 100,
  CostFunction cost_function = CostFunction::PRIMAL_VOLUME);
