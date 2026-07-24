#include "native_hull.h"

#include <array>
#include <utility>

extern "C" {
#include <libqhull_r/qhull_ra.h>
}

namespace nat
{

Mesh convex_hull_3(const std::vector<Point3> & pts)
{
  Mesh M;
  const int n = (int)pts.size();
  if(n < 4) return M;   // no 3D hull

  std::vector<coordT> coords(3 * n);
  for(int i = 0; i < n; ++i)
  {
    coords[3 * i + 0] = pts[i].x();
    coords[3 * i + 1] = pts[i].y();
    coords[3 * i + 2] = pts[i].z();
  }

  qhT qh_qh;
  qhT * qh = &qh_qh;
  qh_zero(qh, stderr);

  // "Qt": triangulate all output facets so every facet is a triangle.
  char flags[] = "qhull Qt";
  const int exitcode =
      qh_new_qhull(qh, 3, n, coords.data(), False, flags, nullptr, stderr);

  std::vector<Point3> hv;                 // compact hull vertices
  std::vector<int> old2new(n, -1);
  std::vector<std::array<int, 3>> tris;

  if(!exitcode)
  {
    facetT * facet;
    vertexT *vertex, **vertexp;
    FORALLfacets
    {
      if(facet->upperdelaunay) continue;

      // Ordered triangle vertices of this (triangulated) facet.
      setT * vs = qh_facet3vertex(qh, facet);
      int ids[3];
      int c = 0;
      FOREACHvertex_(vs)
      {
        if(c < 3) ids[c] = qh_pointid(qh, vertex->point);
        ++c;
      }
      qh_settempfree(qh, &vs);
      if(c != 3) continue;   // skip anything not a triangle

      int m[3];
      for(int k = 0; k < 3; ++k)
      {
        const int pid = ids[k];
        if(old2new[pid] < 0) { old2new[pid] = (int)hv.size(); hv.push_back(pts[pid]); }
        m[k] = old2new[pid];
      }

      // Orient the triangle so its normal agrees with qhull's outward normal.
      const Vector3 nrm = cross_product(pts[ids[1]] - pts[ids[0]],
                                        pts[ids[2]] - pts[ids[0]]);
      const double dp = nrm.x() * facet->normal[0] +
                        nrm.y() * facet->normal[1] +
                        nrm.z() * facet->normal[2];
      if(dp < 0) std::swap(m[1], m[2]);
      tris.push_back({m[0], m[1], m[2]});
    }
  }

  qh_freeqhull(qh, qh_ALL);
  int curlong, totlong;
  qh_memfreeshort(qh, &curlong, &totlong);

  M.build(hv, tris);
  return M;
}

}  // namespace nat
