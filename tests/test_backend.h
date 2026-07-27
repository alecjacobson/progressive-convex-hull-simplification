#pragma once
// Central place for the mesh/kernel types the tests build on. Today this aliases
// the CGAL backend directly; when the backend seam lands (Phase A) this becomes
// `mesh::Polyhedron` / `geom::...` and the same tests validate both backends.
#include "convex_hull_simplification.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

using TestPolyhedron = ConvexHullSimplification::IPolyhedron;

// --- small geometry helpers used by several tests -------------------------

// Regular unit cube [-1,1]^3 corners.
inline Eigen::MatrixXd cube_vertices()
{
  Eigen::MatrixXd V(8, 3);
  int r = 0;
  for(int x : {-1, 1})
    for(int y : {-1, 1})
      for(int z : {-1, 1})
        V.row(r++) << x, y, z;
  return V;
}

// Regular tetrahedron inscribed in the unit cube.
inline Eigen::MatrixXd regular_tetrahedron_vertices()
{
  Eigen::MatrixXd V(4, 3);
  V.row(0) <<  1,  1,  1;
  V.row(1) <<  1, -1, -1;
  V.row(2) << -1,  1, -1;
  V.row(3) << -1, -1,  1;
  return V;
}
