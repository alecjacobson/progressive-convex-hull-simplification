#include "convex_hull_simplification.h"
#include "dual_hull.h"

#include <igl/get_seconds.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/copyleft/cgal/assign.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

// ---------------------------------------------------------------------------
// Constructor: build primal hull → dual hull → initialize queue
// ---------------------------------------------------------------------------

ConvexHullSimplification::ConvexHullSimplification(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int max_degree_for_flips,
  CostFunction cost_function)
  : max_degree_for_flips_(max_degree_for_flips)
  , cost_function_(cost_function)
  , stats_({0, 0, 0, 0})
{
  using EPolyhedron = mesh::Polyhedron;

  // --- Primal convex hull ---
  {
    const double t0 = igl::get_seconds();
    auto primal_points = point_list<EK>(V);
    EPolyhedron primal;
    //printf("primal points: %d\n", (int)primal_points.size());
    mesh::convex_hull_3(primal_points.begin(), primal_points.end(), primal);
    //printf("primal facets: %d\n", (int)primal.size_of_facets());
    stats_.t_primal_hull = igl::get_seconds() - t0;

    // --- Dual hull (Chebyshev center + dual points + dedup + convex_hull_3) ---
    const double t1 = igl::get_seconds();
    const auto bb = mesh::bounding_box(primal_points.begin(), primal_points.end());
    const double bbd_sqr_len = geom::squared_distance(bb.min(), bb.max());
    const double primal_squared_area_tol = 1e-15 * bbd_sqr_len * bbd_sqr_len;
    //printf("primal squared area tol: %g\n", primal_squared_area_tol);
    x0_exact_ = polyhedron_chebyshev_center(primal, primal_squared_area_tol);

    auto dual_points_exact = dual_points_list(primal, primal_squared_area_tol, x0_exact_);
    //printf("Exact dual points: %d\n", (int)dual_points_exact.size());
    auto dV_exact = point_matrix<EK>(dual_points_exact);

    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> dV;
    igl::copyleft::cgal::assign(dV_exact, true, dV);
    {
      Eigen::VectorXi _1, _2;
      igl::remove_duplicate_vertices(
        Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor>(dV),
        1e-14, dV, _1, _2);
    }
    //printf("Deduplicated dual points: %d\n", (int)dV.rows());

    {
      auto dpts = point_list<IK>(dV);
      mesh::convex_hull_3(dpts.begin(), dpts.end(), dual_);
    }
    stats_.t_dual_hull = igl::get_seconds() - t1;
  }
  //printf("dual vertices: %d\n", (int)dual_.size_of_vertices());

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
    pop_ids_.reserve(nv);
    for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
    {
      if(v == nullptr)
      {
        assert(false);
        continue;
      }
      if(v->halfedge() == nullptr)
      {
        // 
#ifndef NDEBUG
        printf("Warning: vertex %d has null halfedge; skipping\n", v->id());
        printf("  %g,%g,%g\n", v->point().x(), v->point().y(), v->point().z());
#endif
        continue;
      }
      auto [cost, record] = measure_vertex_erasure(dual_, v, max_degree_for_flips_, cost_function_);
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
  //printf("\n\n\n\n\n\n\n--------------------------------------\n");
  //printf("mean_width (before): %g  (dual vertices: %d)\n", mean_width(), num_dual_vertices());
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
  pop_ids_.push_back(id);

//#warning "debugging"
  {
    //printf("popping vertex %d with cost %g\n", id, cost);
    auto [debug_cost, debug_record] = measure_vertex_erasure(dual_, v, max_degree_for_flips_, cost_function_);
    //printf("  verified cost: %g\n", debug_cost);
  }


  erase_vertex_and_clip_ears(dual_, v, record);
  //printf("mean_width (after): %g  (dual vertices: %d)\n", mean_width(), num_dual_vertices());

  assert(dual_.is_pure_triangle());
  assert(dual_.is_closed());

  for(const auto & n : neighbors)
  {
    const int nid = n->id();
    auto [ncost, nrecord] = measure_vertex_erasure(dual_, n, max_degree_for_flips_, cost_function_);
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
// dual_vertex_ids(): original vertex ids in dual iteration order
// ---------------------------------------------------------------------------

Eigen::VectorXi ConvexHullSimplification::dual_vertex_ids() const
{
  // Iterate in the same order as dual_to_primal_mesh so that result[i]
  // is the id of the dual vertex that produces polygon i in pPI/pPC.
  Eigen::VectorXi ids(num_dual_vertices());
  int i = 0;
  for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
    ids(i++) = static_cast<int>(v->id());
  return ids;
}

// ---------------------------------------------------------------------------
// dual_greedy_coloring(): greedy graph coloring of the dual's 1-skeleton
// ---------------------------------------------------------------------------

Eigen::VectorXi ConvexHullSimplification::dual_greedy_coloring(int k) const
{
  // Output indexed by vertex id; -1 means the vertex has been removed.
  const int N = static_cast<int>(full_records_.size());
  Eigen::VectorXi colors = Eigen::VectorXi::Constant(N, -1);

  // --- Initial greedy assignment: smallest available label in [0, k) ---
  for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
  {
    const int vid = static_cast<int>(v->id());
    std::vector<bool> used(k, false);
    auto h = v->halfedge();
    do {
      const int nid = static_cast<int>(h->opposite()->vertex()->id());
      if(colors(nid) >= 0) used[colors(nid)] = true;
      h = h->next()->opposite();
    } while(h != v->halfedge());
    int c = 0;
    while(c < k-1 && used[c]) ++c;  // cap at k-1; label is always in [0, k)
    colors(vid) = c;
  }

  // --- Rebalancing post-process ---
  // Repeatedly try to move vertices from over-populated to under-populated
  // labels. A move is accepted only when the count gap is >= 2, which
  // strictly decreases sum(cnt[c]^2) and guarantees termination.
  std::vector<int> cnt(k, 0);
  for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
    ++cnt[colors(static_cast<int>(v->id()))];

  bool changed = true;
  while(changed)
  {
    changed = false;
    for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
    {
      const int vid = static_cast<int>(v->id());
      const int cur = colors(vid);

      std::vector<bool> forbidden(k, false);
      auto h = v->halfedge();
      do {
        const int nid = static_cast<int>(h->opposite()->vertex()->id());
        if(colors(nid) >= 0) forbidden[colors(nid)] = true;
        h = h->next()->opposite();
      } while(h != v->halfedge());

      // Pick the available label with lowest count that is at least 2 below
      // cnt[cur], so the sum-of-squares potential strictly decreases.
      int best = -1;
      for(int c = 0; c < k; ++c)
      {
        if(forbidden[c]) continue;
        if(cnt[c] + 1 < cnt[cur] && (best == -1 || cnt[c] < cnt[best]))
          best = c;
      }
      if(best != -1)
      {
        --cnt[cur];
        ++cnt[best];
        colors(vid) = best;
        changed = true;
      }
    }
  }

  // check all edges are properly colored
  for(auto e = dual_.edges_begin(); e != dual_.edges_end(); ++e)
  {
    auto v1 = e->vertex();
    auto v2 = e->opposite()->vertex();
    if(colors(static_cast<int>(v1->id())) == colors(static_cast<int>(v2->id())))
    {
      throw std::runtime_error("Coloring failed: adjacent vertices " +
        std::to_string(v1->id()) + " and " + std::to_string(v2->id()) +
        " share color " + std::to_string(colors(static_cast<int>(v1->id()))));
    }
  }

  return colors;
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
// get_dual_mesh(): extract dual as a triangle (V, F) with remapped indices
// ---------------------------------------------------------------------------

std::pair<Eigen::MatrixXd, Eigen::MatrixXi>
ConvexHullSimplification::get_dual_mesh() const
{
  const int nv = num_dual_vertices();
  Eigen::MatrixXd V(nv, 3);
  std::vector<int> id_to_idx(full_records_.size(), -1);
  {
    int vi = 0;
    for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
    {
      const int id = static_cast<int>(v->id());
      id_to_idx[id] = vi;
      V.row(vi) << v->point().x(), v->point().y(), v->point().z();
      vi++;
    }
  }

  assert(dual_.is_pure_triangle());
  const int nf = static_cast<int>(dual_.size_of_facets());
  Eigen::MatrixXi F(nf, 3);
  {
    int fi = 0;
    for(auto f = dual_.facets_begin(); f != dual_.facets_end(); ++f)
    {
      auto h = f->halfedge();
      F.row(fi++) <<
        id_to_idx[h->vertex()->id()],
        id_to_idx[h->next()->vertex()->id()],
        id_to_idx[h->next()->next()->vertex()->id()];
    }
  }
  return {V, F};
}

Eigen::VectorXd ConvexHullSimplification::popped_dual_vertex_costs() const
{
  Eigen::VectorXd costs(full_records_.size());
  for(size_t i = 0; i < full_records_.size(); ++i)
  {
    costs(i) = full_records_[i].cost;
  }
  for(auto v = dual_.vertices_begin(); v != dual_.vertices_end(); ++v)
  {
    const int vid = static_cast<int>(v->id());
    costs(vid) = std::numeric_limits<Scalar>::quiet_NaN();
  }
  return costs;
}

Eigen::VectorXi ConvexHullSimplification::popped_dual_vertex_ids() const
{
  return Eigen::VectorXi(Eigen::Map<const Eigen::VectorXi>(pop_ids_.data(), pop_ids_.size()));
}


// ---------------------------------------------------------------------------
// stats()
// ---------------------------------------------------------------------------

ConvexHullSimplification::Stats ConvexHullSimplification::stats() const
{
  return stats_;
}

// ---------------------------------------------------------------------------
// mean_width(): (1/4π) Σ ℓ_e θ_e over all primal edges
// ---------------------------------------------------------------------------

double ConvexHullSimplification::mean_width() const
{
  // Chebyshev center (EK == IK here, both are inexact constructions kernel)
  const double x0x = geom::to_double(x0_exact_.x());
  const double x0y = geom::to_double(x0_exact_.y());
  const double x0z = geom::to_double(x0_exact_.z());

  // For a dual triangular face, return the corresponding primal vertex position
  // via the inverse polarity map: p = x0 + n_face / (n_face · centroid).
  const auto primal_vertex = [&](auto f) -> Eigen::Vector3d
  {
    auto h = f->halfedge();
    const auto & A = h->vertex()->point();
    const auto & B = h->next()->vertex()->point();
    const auto & C = h->next()->next()->vertex()->point();
    const auto n  = geom::cross_product(B - A, C - A);
    const auto bc = ((A - geom::ORIGIN) + (B - geom::ORIGIN) + (C - geom::ORIGIN)) / IK::FT(3);
    const double beta = geom::to_double(n * bc);
    return {x0x + geom::to_double(n.x()) / beta,
            x0y + geom::to_double(n.y()) / beta,
            x0z + geom::to_double(n.z()) / beta};
  };

  double sum = 0.0;
  for(auto e = dual_.edges_begin(); e != dual_.edges_end(); ++e)
  {
    // Dual vertex position ∝ outward primal face normal (dual point = n/|denom|, denom < 0)
    const auto & p1 = e->vertex()->point();
    const auto & p2 = e->opposite()->vertex()->point();
    // These don't need to be normalized if we're using the atan2 version
    const Eigen::Vector3d n1 = Eigen::Vector3d(
      geom::to_double(p1.x()), geom::to_double(p1.y()), geom::to_double(p1.z()));
    const Eigen::Vector3d n2 = Eigen::Vector3d(
      geom::to_double(p2.x()), geom::to_double(p2.y()), geom::to_double(p2.z()));

    // Dual faces give the two endpoints of the primal edge
    const Eigen::Vector3d q0 = primal_vertex(e->facet());
    const Eigen::Vector3d q1 = primal_vertex(e->opposite()->facet());

    const Eigen::Vector3d edge_vec = q1 - q0;
    const double len   = (q1 - q0).stableNorm();
    if(len == 0.0)
    {
      continue;
    }
    const Eigen::Vector3d edge_dir = edge_vec / len;
    //const double theta = std::acos(std::clamp(n1.dot(n2), -1.0, 1.0));
    const double theta = std::atan2(
      edge_dir.dot(n1.cross(n2)),
      n1.dot(n2));
    const auto contrib = len * theta;
    //printf("%d,%d  %g\n",
    //    std::min(static_cast<int>(e->vertex()->id()), static_cast<int>(e->opposite()->vertex()->id())),
    //    std::max(static_cast<int>(e->vertex()->id()), static_cast<int>(e->opposite()->vertex()->id())),
    //    contrib);
    sum += contrib;
  }

  return sum / (4.0 * std::acos(-1.0));
}

// ---------------------------------------------------------------------------
// Free-function wrapper
// ---------------------------------------------------------------------------

std::tuple<Eigen::MatrixXd, Eigen::VectorXi, Eigen::VectorXi>
simplify_convex_hull(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int target_num_dual_vertices,
  int max_degree_for_flips,
  CostFunction cost_function)
{
  ConvexHullSimplification chs(V, F, max_degree_for_flips, cost_function);
  chs.simplify_to(target_num_dual_vertices);
  return chs.get_primal_mesh();
}
