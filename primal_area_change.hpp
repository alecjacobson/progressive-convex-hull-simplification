#pragma once
#include <Eigen/Core>
#include <tuple>

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, typename Polyhedron::Traits::FT, bool> primal_area_change(
  const Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0);
