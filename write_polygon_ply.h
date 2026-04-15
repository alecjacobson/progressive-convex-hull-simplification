#pragma once
// Write a polygon mesh (pV, pPI, pPC) as a binary little-endian PLY file.
// Mirrors the minimal igl::writePLY(filename, V, F) API but accepts the
// libigl polygon-soup format instead of a triangle matrix.
//
// pV  — #vertices x 3 vertex positions
// pPI — flat polygon vertex index list
// pPC — polygon cumulative counts, size = num_polygons + 1
//
// Returns true on success.
//
// Note: tinyply does not support variable-length list writing, so this is
// written directly against the PLY spec.

#include <Eigen/Core>
#include <fstream>
#include <string>
#include <cstdint>

inline bool write_polygon_ply(
  const std::string        & filename,
  const Eigen::MatrixXd    & pV,
  const Eigen::VectorXi    & pPI,
  const Eigen::VectorXi    & pPC)
{
  std::ofstream out(filename, std::ios::binary);
  if(!out) return false;

  const int nv = pV.rows();
  const int nf = pPC.size() - 1;

  // ASCII header
  out << "ply\n"
      << "format binary_little_endian 1.0\n"
      << "element vertex " << nv << "\n"
      << "property double x\n"
      << "property double y\n"
      << "property double z\n"
      << "element face " << nf << "\n"
      << "property list uchar int vertex_indices\n"
      << "end_header\n";

  // Vertices — 3 doubles per vertex
  for(int i = 0; i < nv; i++)
  {
    double xyz[3] = { pV(i,0), pV(i,1), pV(i,2) };
    out.write(reinterpret_cast<const char*>(xyz), sizeof(xyz));
  }

  // Faces — [uint8 count, int32 v0, ..., int32 v_{count-1}]
  for(int p = 0; p < nf; p++)
  {
    const int start = pPC(p);
    const uint8_t cnt = static_cast<uint8_t>(pPC(p+1) - pPC(p));
    out.write(reinterpret_cast<const char*>(&cnt), 1);
    for(int k = 0; k < cnt; k++)
    {
      int32_t idx = static_cast<int32_t>(pPI(start + k));
      out.write(reinterpret_cast<const char*>(&idx), 4);
    }
  }

  return out.good();
}
