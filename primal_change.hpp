#pragma once
#include "CostFunction.hpp"
#include "vertex_erasure.h"
#include <tuple>

/// poly  a copy of the dual one-ring after the vertex that used to correspond
/// to p0 has been removed and the hole has been re-triangulated convexly.
template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, bool> primal_change(
  Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0,
  const CostFunction cost_function = CostFunction::PRIMAL_VOLUME);
