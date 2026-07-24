#pragma once
// -----------------------------------------------------------------------------
// Backend selection: geometry kernel (points/vectors + predicates + constructions).
//
// Second of the two seam headers (the other is mesh_backend.h). Every geometry
// primitive the algorithm uses is named through `namespace geom`, so a backend
// switch only changes what these names resolve to:
//
//   * PCHS_BACKEND_CGAL (default) — forwards to CGAL's exact/inexact predicates
//     and constructions.
//   * PCHS_BACKEND_NATIVE (Phase B/C) — native double point/vector types with
//     Shewchuk's robust orient3d and plain double arithmetic for the rest.
//
// The ONLY predicate that needs robustness is orient3d (drives every convexity
// decision); the rest are ordinary double arithmetic.
// -----------------------------------------------------------------------------

#include "mesh_backend.h"   // pulls PCHS_BACKEND_* selection

#if defined(PCHS_BACKEND_CGAL)

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>            // to_double, sqrt
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/squared_distance_3.h>

namespace geom
{
  // Predicates / constructions.
  using CGAL::orientation;       // orient3d: Sign of (a,b,c,d)
  using CGAL::cross_product;
  using CGAL::scalar_product;
  using CGAL::centroid;
  using CGAL::squared_area;
  using CGAL::squared_distance;
  using CGAL::to_double;
  using CGAL::sqrt;

  // Constants.
  using CGAL::ORIGIN;
  using CGAL::NULL_VECTOR;

  // Sign / oriented-side enumerators.
  using CGAL::NEGATIVE;
  using CGAL::ZERO;
  using CGAL::POSITIVE;
  using CGAL::ON_NEGATIVE_SIDE;
  using CGAL::ON_ORIENTED_BOUNDARY;
  using CGAL::ON_POSITIVE_SIDE;
}

#elif defined(PCHS_BACKEND_NATIVE)

#include "native_geom.h"

namespace geom
{
  using nat::orientation;
  using nat::cross_product;
  using nat::scalar_product;
  using nat::centroid;
  using nat::squared_area;
  using nat::squared_distance;
  using nat::to_double;
  using nat::sqrt;

  using nat::ORIGIN;
  using nat::NULL_VECTOR;

  using nat::NEGATIVE;
  using nat::ZERO;
  using nat::POSITIVE;
  using nat::ON_NEGATIVE_SIDE;
  using nat::ON_ORIENTED_BOUNDARY;
  using nat::ON_POSITIVE_SIDE;
}

#endif
