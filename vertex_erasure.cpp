#include "vertex_erasure.h"
#include "primal_change.hpp"
#include "dual_hull.h"
#include "polyhedron_utils.h"

#include <cassert>
#include <stdexcept>


// Ear-clip one_ring_copy (a triangulated disk with open boundary) starting
// from h0. Records the clipping order as a path (walks between ears) and the
// starting vertex id so it can be replayed with erase_vertex_and_clip_ears.
template <class Polyhedron>
std::pair<size_t, std::vector<int>>
clip_ears(Polyhedron & one_ring_copy, typename Polyhedron::Halfedge_handle h0)
{
  const auto walk = [](const auto h)->auto
  {
    assert(h->opposite()->is_border());
    return h->opposite()->prev()->opposite();
  };

#ifndef NDEBUG
  {
    int count = 0;
    auto h = h0;
    while(true)
    {
      h = walk(h);
      count++;
      assert(count < 1000);
      if(h == h0) break;
    }
    assert(count == (int)one_ring_copy.size_of_vertices());
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
    if(w == n)
    {
      // Ear found
      path.push_back(walks_since_last_ear);
      walks_since_last_ear = 0;
      auto q = h->next()->next()->opposite();
      one_ring_copy.erase_facet(h);
      assert(q->opposite()->is_border());
      h = q;
    }
    else
    {
      h = w;
      walks_since_last_ear++;
    }
  }
  return {start_vertex_id, path};
}

// Apply a previously recorded ear-clipping path to the actual dual polyhedron.
// Erases v and re-triangulates the resulting hole using the recorded path.
template <class Polyhedron>
void erase_vertex_and_clip_ears(
  Polyhedron & dual,
  typename Polyhedron::Vertex_handle & v,
  const size_t start_vertex_id,
  const std::vector<int> & path)
{
  auto h0 = dual.erase_center_vertex(v->halfedge());
  while(h0->vertex()->id() != start_vertex_id)
    h0 = h0->next();

  const auto walk = [](const auto h)->auto { return h->next(); };

  auto h = h0;
  for(int i = 0; i < (int)path.size(); i++)
  {
    const int steps = path[i];
    for(int j = 0; j < steps; j++)
      h = walk(h);
    h = dual.split_facet(h->prev(), h->next());
  }
}

// Collect all neighbor vertex handles of v in one-ring order.
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

// Compute the cost (primal volume added) of erasing dual vertex v, and record
// the ear-clipping path needed to re-triangulate its hole.
template <class Polyhedron>
std::pair<typename Polyhedron::Traits::FT, Record>
measure_vertex_erasure(
  const Polyhedron & dual,
  const typename Polyhedron::Vertex_handle & v,
  const int max_degree_for_flips,
  const CostFunction cost_function)
{
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

  if((int)v->degree() <= max_degree_for_flips)
    one_ring_triangulation_convex_via_flips(dual,v,one_ring_copy,h0);
  else
    one_ring_triangulation_convex_via_convex_hull(dual,v,one_ring_copy,h0);

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
      throw std::runtime_error("Removed point lies on negative side of a new face.");
  }
#endif
//      if(v->id() == 10)
//      {
//        printf("------------------------------------------------------\n");
//
//        print_VF(dual,"d");
//        print_VF(one_ring_copy,"o");
//        using Point = typename Polyhedron::Traits::Point_3;
//        const auto print_primal = [](const auto & dual, const std::string prefix)
//        {
//          Point x0(0.,0.,0.);
//          {
//            Polyhedron dual_cpy = dual;
//            {
//              int fid = 0;
//              for(auto f = dual_cpy.facets_begin(); f != dual_cpy.facets_end(); ++f)
//                f->id() = fid++;
//            }
//            auto [pV, pPI, pPC] = dual_to_primal_mesh(dual_cpy, x0);
//            std::cout<<igl::matlab_format(pV,prefix+"V")<<std::endl;
//            std::cout<<igl::matlab_format_index(pPI.transpose().eval(),prefix+"PI")<<std::endl;
//            std::cout<<igl::matlab_format(pPC.transpose().eval(),prefix+"PC")<<std::endl;
//          }
//        };
//        const auto print_primal_points = [](const auto & dual, const std::string prefix)
//        {
//          Point x0(0.,0.,0.);
//          {
//            Polyhedron dual_cpy = dual;
//            {
//              int fid = 0;
//              for(auto f = dual_cpy.facets_begin(); f != dual_cpy.facets_end(); ++f)
//                f->id() = fid++;
//            }
//            auto pV = dual_to_primal_points(dual_cpy, x0);
//            std::cout<<igl::matlab_format(pV,prefix+"V")<<std::endl;
//          }
//        };
//        print_primal(dual,"p");
//        print_primal_points(one_ring_copy,"q");
//        printf(R"(id = %d+1;
//rF = [dF(any(dF==id,2),:);fliplr(oF)];
//rN = normals(dV,rF);
//rBC = barycenter(dV,rF);
//rbeta = sum(rN.*rBC,2);
//qV = rN./rbeta;
//[~,qPI,qPC] = dual(dV,rF)
//[qF,qJ] = polygons_to_triangles(qPI,qPC);
//)",v->id());
//      }

  

  auto [cost,contains_origin] = primal_change(
    one_ring_copy,
    p,
    cost_function);
  //printf("cost(%d): %g\n",v->id(),cost);
  //    if(v->id() == 10)
  //    {
  //      exit(1);
  //    }

  Record record;
  std::tie(record.start_vertex_id, record.path) = clip_ears(one_ring_copy, h0);

  return {cost, record};
}


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
template std::vector<CGAL::internal::In_place_list_iterator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>, std::allocator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>>>, std::allocator<CGAL::internal::In_place_list_iterator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>, std::allocator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>>>>> collect_neighbors<CGAL::internal::In_place_list_iterator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>, std::allocator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>>>>(CGAL::internal::In_place_list_iterator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>, std::allocator<CGAL::HalfedgeDS_in_place_list_vertex<CGAL::I_Polyhedron_vertex<CGAL::HalfedgeDS_vertex_max_base_with_id<CGAL::HalfedgeDS_list_types<CGAL::Epick, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int>>, CGAL::Point_3<CGAL::Epick>, unsigned long>>>>> const&);
template std::pair<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Traits::FT, Record> measure_vertex_erasure<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>> const&, CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Vertex_handle const&, int, CostFunction);
template void erase_vertex_and_clip_ears<CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>>(CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>&, CGAL::Polyhedron_3<CGAL::Epick, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int>>::Vertex_handle&, unsigned long, std::vector<int, std::allocator<int>> const&);

