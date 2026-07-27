// Robust orient3d for the native backend, via Shewchuk's adaptive exact
// predicates (danshapero/predicates). Compiled into pchs_predicates_wrap, which
// defines PCHS_ROBUST_PREDICATES so native_geom.h routes nat::orient3d here.
#include "native_geom.h"

#include <predicates.h>   // Shewchuk orient3d (C, extern "C" in the header)

namespace nat
{
namespace detail
{

Sign robust_orient3d(const Point3 & a, const Point3 & b,
                     const Point3 & c, const Point3 & d)
{
  double pa[3] = {a.x(), a.y(), a.z()};
  double pb[3] = {b.x(), b.y(), b.z()};
  double pc[3] = {c.x(), c.y(), c.z()};
  double pd[3] = {d.x(), d.y(), d.z()};
  // Shewchuk's orient3d sign is the negation of CGAL::orientation's convention
  // (verified: it returns -1 for the canonical POSITIVE tetrahedron), so flip.
  const double r = predicates::orient3d(pa, pb, pc, pd);   // Shewchuk
  if(r < 0) return POSITIVE;
  if(r > 0) return NEGATIVE;
  return ZERO;
}

}  // namespace detail
}  // namespace nat
