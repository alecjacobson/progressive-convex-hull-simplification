#pragma once
#include <Eigen/Core>
#include <tuple>

/// @param[in]  V  #V by 3 list of dual vertex positions
template 
  <typename DerivedV,
   typename DerivedPI,
   typename DerivedPC>
std::tuple<typename DerivedV::Scalar, typename DerivedV::Scalar, bool> primal_volume_subtended(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedPI> & PI,
  const Eigen::MatrixBase<DerivedPC> & PC,
  const int i0);

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, typename Polyhedron::Traits::FT, bool> primal_volume_subtended(
  const Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0);

template 
  <typename DerivedV>
std::tuple<typename DerivedV::Scalar, typename DerivedV::Scalar, bool> primal_volume_subtended(
  const Eigen::MatrixBase<DerivedV> & V);

template 
  <typename Polyhedron>
std::tuple<typename Polyhedron::Traits::FT, bool> primal_volume_subtended_explicit(
  const Polyhedron & poly,
  const typename Polyhedron::Point_3 & p0);
