#pragma once
#include "vertex_erasure.h"

#include <igl/min_heap.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

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
  using EK = CGAL::Exact_predicates_exact_constructions_kernel;
  using IK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using IPolyhedron = CGAL::Polyhedron_3<IK, CGAL::Polyhedron_items_with_id_3>;
  using Scalar = IK::FT;

  // Build the primal convex hull, compute the Chebyshev center, construct the
  // dual, and initialize the priority queue with costs for every dual vertex.
  ConvexHullSimplification(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    int max_degree_for_flips = 100);

  // Drain stale queue entries, then pop and apply the cheapest non-stale
  // vertex erasure. Updates neighbor costs. Returns false if the queue is
  // exhausted or every remaining vertex has infinite cost.
  bool step();

  // Repeatedly call step() until num_dual_vertices() == target or step()
  // returns false. No-op if already at or below target.
  void simplify_to(int target_num_dual_vertices);

  // Number of dual vertices currently remaining (= number of halfspaces).
  int num_dual_vertices() const;

  // Reconstruct the primal polygon mesh from the current dual state.
  // pV  — #dual_facets x 3 primal vertex positions (one per dual face)
  // pPI — flat polygon vertex index list
  // pPC — cumulative counts (size = num_dual_vertices() + 1)
  std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi> get_primal_mesh();

  struct Stats
  {
    double t_primal_hull;    // convex_hull_3 on input points
    double t_dual_hull;      // Chebyshev center + dual points + dedup + convex_hull_3
    double t_queue_init;     // initial measure_vertex_erasure for all dual vertices
    double t_last_simplify;  // most recent simplify_to() call (0 if never called)
  };
  Stats stats() const;

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
  igl::min_heap<std::tuple<Scalar,int,int>> Q_;
  int max_degree_for_flips_;
  Stats stats_;
};

// Convenience wrapper: construct, simplify to target, return polygon mesh.
std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi>
simplify_convex_hull(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int target_num_dual_vertices,
  int max_degree_for_flips = 100);
