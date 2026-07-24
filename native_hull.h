#pragma once
// Native convex hull (Backend B): qhull -> nat::Mesh. Std-only interface; the
// qhull dependency is confined to native_hull.cpp.
#include "native_mesh.h"

#include <vector>

namespace nat
{
// Convex hull of pts as a closed, outward-oriented triangle mesh. Only hull
// vertices appear (interior points are dropped); facets are triangulated (Qt).
Mesh convex_hull_3(const std::vector<Point3> & pts);
}
