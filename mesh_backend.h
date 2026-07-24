#pragma once
// -----------------------------------------------------------------------------
// Backend selection: mesh datastructure + global hull/bounding-box algorithms.
//
// This is one of the two seam headers (the other is geometry.h). Everything the
// algorithm needs from the half-edge datastructure and the global geometry
// algorithms is named through `namespace mesh`. Selecting a backend is a matter
// of switching what these names alias:
//
//   * PCHS_BACKEND_CGAL (default) — CGAL::Polyhedron_3 + CGAL::convex_hull_3.
//   * PCHS_BACKEND_NATIVE (Phase B/C) — std-only half-edge mesh + qhull.
//
// The mesh type deliberately mirrors the CGAL Polyhedron_3 member API
// (Halfedge_handle/->next()/->opposite()/Euler ops/...), so the templated
// algorithm code is instantiated on `mesh::Polyhedron` with no other changes.
// -----------------------------------------------------------------------------

#if !defined(PCHS_BACKEND_NATIVE)
#  define PCHS_BACKEND_CGAL 1
#endif

#if defined(PCHS_BACKEND_CGAL)

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/bounding_box.h>

namespace mesh
{
  // Single inexact (double) kernel throughout — the exact-kernel path was
  // already disabled in the original code.
  using Kernel     = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point3     = Kernel::Point_3;
  using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;

  // Convex hull of a point range into a half-edge mesh.
  template <class InputIt, class Poly>
  inline void convex_hull_3(InputIt begin, InputIt end, Poly & out)
  {
    CGAL::convex_hull_3(begin, end, out);
  }

  // Axis-aligned bounding box of a point range (as CGAL::Iso_cuboid_3).
  template <class InputIt>
  inline auto bounding_box(InputIt begin, InputIt end)
  {
    return CGAL::bounding_box(begin, end);
  }
}

#elif defined(PCHS_BACKEND_NATIVE)
#  error "Native backend (Phase B/C) not wired up yet."
#endif
