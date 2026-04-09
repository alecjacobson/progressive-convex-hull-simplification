#include "convex_hull_simplification.h"
#include "dual_hull.h"

#include <igl/get_seconds.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/copyleft/cgal/assign.h>

#include <CGAL/convex_hull_3.h>

#include <cassert>
#include <limits>

// ---------------------------------------------------------------------------
// Constructor: build primal hull → dual hull → initialize queue
// ---------------------------------------------------------------------------

ConvexHullSimplification::ConvexHullSimplification(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int max_degree_for_flips)
  : max_degree_for_flips_(max_degree_for_flips)
  , stats_({0, 0, 0, 0})
{
  using EPolyhedron = CGAL::Polyhedron_3<EK, CGAL::Polyhedron_items_with_id_3>;

  // --- Primal convex hull (exact kernel) ---
  {
    const double t0 = igl::get_seconds();
    auto primal_points = point_list<EK>(V);
    EPolyhedron primal;
    CGAL::convex_hull_3(primal_points.begin(), primal_points.end(), primal);
    stats_.t_primal_hull = igl::get_seconds() - t0;

    // --- Dual hull (Chebyshev center + dual points + dedup + convex_hull_3) ---
    const double t1 = igl::get_seconds();
    const double primal_squared_area_tol = 1e-15;
    x0_exact_ = polyhedron_chebyshev_center(primal, primal_squared_area_tol);

    auto dual_points_exact = dual_points_list(primal, primal_squared_area_tol, x0_exact_);
    auto dV_exact = point_matrix<EK>(dual_points_exact);

    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> dV;
    igl::copyleft::cgal::assign(dV_exact, true, dV);
    {
      Eigen::VectorXi _1, _2;
      igl::remove_duplicate_vertices(
        Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor>(dV),
        1e-14, dV, _1, _2);
    }

    {
      auto dpts = point_list<IK>(dV);
      CGAL::convex_hull_3(dpts.begin(), dpts.end(), dual_);
    }
    stats_.t_dual_hull = igl::get_seconds() - t1;
  }

  // Assign contiguous vertex ids
  {
    int vid = 0;
    for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
      v->id() = vid++;
  }

#ifndef NDEBUG
  confirm_all_edges_are_convex(dual_);
#endif

  assert(dual_.is_pure_triangle());
  assert(dual_.is_closed());

  // --- Initial queue fill ---
  {
    const double t0 = igl::get_seconds();
    const int nv = dual_.size_of_vertices();
    full_records_.resize(nv);
    for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
    {
      auto [cost, record] = measure_vertex_erasure(dual_, v, max_degree_for_flips_);
      const int vid = v->id();
      full_records_[vid] = {cost, 0, v, record};
      Q_.push({cost, vid, 0});
    }
    stats_.t_queue_init = igl::get_seconds() - t0;
  }
}

// ---------------------------------------------------------------------------
// step(): drain stale entries, apply cheapest non-stale erasure
// ---------------------------------------------------------------------------

bool ConvexHullSimplification::step()
{
  // Drain stale entries
  while(!Q_.empty())
  {
    auto [cost, id, visit] = Q_.top();
    if(visit == full_records_[id].visit) break;
    Q_.pop();
  }

  if(Q_.empty()) return false;

  auto [cost, id, visit] = Q_.top();
  Q_.pop();

  if(cost == std::numeric_limits<Scalar>::infinity()) return false;

  auto v = full_records_[id].vertex;
  const auto & record = full_records_[id].record;
  const auto neighbors = collect_neighbors(v);

  erase_vertex_and_clip_ears(dual_, v, record.start_vertex_id, record.path);

  assert(dual_.is_pure_triangle());
  assert(dual_.is_closed());

  for(const auto & n : neighbors)
  {
    const int nid = n->id();
    auto [ncost, nrecord] = measure_vertex_erasure(dual_, n, max_degree_for_flips_);
    full_records_[nid].cost   = ncost;
    full_records_[nid].record = nrecord;
    full_records_[nid].visit++;
    Q_.push({ncost, nid, full_records_[nid].visit});
  }

  return true;
}

// ---------------------------------------------------------------------------
// simplify_to(): drive step() until target or exhausted
// ---------------------------------------------------------------------------

void ConvexHullSimplification::simplify_to(int target_num_dual_vertices)
{
  const double t0 = igl::get_seconds();
  while(num_dual_vertices() > target_num_dual_vertices)
  {
    if(!step()) break;
  }
  stats_.t_last_simplify = igl::get_seconds() - t0;
}

// ---------------------------------------------------------------------------
// num_dual_vertices()
// ---------------------------------------------------------------------------

int ConvexHullSimplification::num_dual_vertices() const
{
  return static_cast<int>(dual_.size_of_vertices());
}

// ---------------------------------------------------------------------------
// get_primal_mesh(): assign face ids, convert dual → primal polygon mesh
// ---------------------------------------------------------------------------

std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi>
ConvexHullSimplification::get_primal_mesh()
{
  // (Re)assign face ids — new faces from split_facet start at 0, so always redo
  {
    int fid = 0;
    for(auto f = dual_.facets_begin(); f != dual_.facets_end(); ++f)
      f->id() = fid++;
  }

  auto [pV_exact, pPI, pPC] = dual_to_primal_mesh(dual_, x0_exact_);

  Eigen::MatrixXd pV;
  igl::copyleft::cgal::assign(pV_exact, pV);

  return {pV, pPI, pPC};
}

// ---------------------------------------------------------------------------
// stats()
// ---------------------------------------------------------------------------

ConvexHullSimplification::Stats ConvexHullSimplification::stats() const
{
  return stats_;
}

// ---------------------------------------------------------------------------
// Free-function wrapper
// ---------------------------------------------------------------------------

std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi>
simplify_convex_hull(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int target_num_dual_vertices,
  int max_degree_for_flips)
{
  ConvexHullSimplification chs(V, F, max_degree_for_flips);
  chs.simplify_to(target_num_dual_vertices);
  return chs.get_primal_mesh();
}
