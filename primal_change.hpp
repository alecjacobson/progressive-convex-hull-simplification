#pragma once
#include "CostFunction.hpp"
#include "vertex_erasure.h"
#include <tuple>

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, bool> primal_change(
  Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0,
  const CostFunction cost_function = CostFunction::PRIMAL_VOLUME);
