#pragma once
// -----------------------------------------------------------------------------
// Native geometry kernel (Backend B) — std-only double point/vector types and
// the geometry primitives the algorithm needs, mirroring the CGAL API so that
// geometry.h's PCHS_BACKEND_NATIVE branch can `using`-import these into `geom`
// and the templated algorithm compiles on the native mesh unchanged.
//
// Point3/Vector3 expose CGAL-style .x()/.y()/.z() accessors and
// (double,double,double) constructors. Everything is plain double arithmetic
// EXCEPT orient3d, the one predicate whose sign drives every convexity
// decision; it currently uses a straight determinant (exact for integer/dyadic
// inputs, which is what the unit tests use). Phase-B hardening replaces its body
// with Shewchuk's adaptive robust orient3d — call sites do not change.
// -----------------------------------------------------------------------------
#include <cmath>

namespace nat
{

// Tag types mirroring CGAL::ORIGIN / CGAL::NULL_VECTOR (defined early so the
// vector/point constructors can accept them).
struct Origin_t {};
struct Null_vector_t {};

class Vector3
{
public:
  Vector3() : x_(0), y_(0), z_(0) {}
  Vector3(double x, double y, double z) : x_(x), y_(y), z_(z) {}
  Vector3(Null_vector_t) : x_(0), y_(0), z_(0) {}
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }
  double squared_length() const { return x_ * x_ + y_ * y_ + z_ * z_; }
private:
  double x_, y_, z_;
};

class Point3
{
public:
  Point3() : x_(0), y_(0), z_(0) {}
  Point3(double x, double y, double z) : x_(x), y_(y), z_(z) {}
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }
  bool operator==(const Point3 & o) const
  { return x_ == o.x_ && y_ == o.y_ && z_ == o.z_; }
  bool operator!=(const Point3 & o) const { return !(*this == o); }
  // Lexicographic order so Point3 can key a std::map (used by the convex-hull
  // one-ring triangulation).
  bool operator<(const Point3 & o) const
  {
    if(x_ != o.x_) return x_ < o.x_;
    if(y_ != o.y_) return y_ < o.y_;
    return z_ < o.z_;
  }
private:
  double x_, y_, z_;
};

// Kernel-like traits bundle (mirrors the CGAL kernel's nested types used by the
// algorithm: Point_3, Vector_3, FT).
struct Kernel
{
  using FT = double;
  using Point_3 = Point3;
  using Vector_3 = Vector3;
};

// Axis-aligned bounding box with CGAL-Iso_cuboid-like min()/max().
struct Bbox
{
  Point3 mn, mx;
  const Point3 & min() const { return mn; }
  const Point3 & max() const { return mx; }
};

inline constexpr Origin_t      ORIGIN{};
inline constexpr Null_vector_t NULL_VECTOR{};

// Sign / oriented-side enumeration mirroring CGAL::Sign (values match).
enum Sign { NEGATIVE = -1, ZERO = 0, POSITIVE = 1 };
inline constexpr Sign ON_NEGATIVE_SIDE     = NEGATIVE;
inline constexpr Sign ON_ORIENTED_BOUNDARY = ZERO;
inline constexpr Sign ON_POSITIVE_SIDE     = POSITIVE;

// --- affine operators (match CGAL's point/vector algebra) -----------------

inline Vector3 operator-(const Point3 & a, const Point3 & b)
{ return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()}; }

inline Vector3 operator-(const Point3 & a, Origin_t)
{ return {a.x(), a.y(), a.z()}; }

inline Point3 operator+(Origin_t, const Vector3 & v)
{ return {v.x(), v.y(), v.z()}; }

inline Point3 operator+(const Point3 & a, const Vector3 & v)
{ return {a.x() + v.x(), a.y() + v.y(), a.z() + v.z()}; }

inline Vector3 operator+(const Vector3 & a, const Vector3 & b)
{ return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()}; }

inline Vector3 operator-(const Vector3 & a, const Vector3 & b)
{ return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()}; }

inline Vector3 operator-(const Vector3 & v)
{ return {-v.x(), -v.y(), -v.z()}; }

inline Vector3 operator*(double s, const Vector3 & v)
{ return {s * v.x(), s * v.y(), s * v.z()}; }

inline Vector3 operator*(const Vector3 & v, double s)
{ return {s * v.x(), s * v.y(), s * v.z()}; }

inline Vector3 operator/(const Vector3 & v, double s)
{ return {v.x() / s, v.y() / s, v.z() / s}; }

// CGAL spells the dot product as operator* on two vectors.
inline double operator*(const Vector3 & a, const Vector3 & b)
{ return a.x() * b.x() + a.y() * b.y() + a.z() * b.z(); }

// --- named primitives (mirror CGAL free functions) ------------------------

inline Vector3 cross_product(const Vector3 & a, const Vector3 & b)
{
  return {a.y() * b.z() - a.z() * b.y(),
          a.z() * b.x() - a.x() * b.z(),
          a.x() * b.y() - a.y() * b.x()};
}

inline double scalar_product(const Vector3 & a, const Vector3 & b)
{ return a * b; }

inline Point3 centroid(const Point3 & a, const Point3 & b, const Point3 & c)
{ return {(a.x() + b.x() + c.x()) / 3.0,
          (a.y() + b.y() + c.y()) / 3.0,
          (a.z() + b.z() + c.z()) / 3.0}; }

inline double squared_area(const Point3 & a, const Point3 & b, const Point3 & c)
{
  const Vector3 n = cross_product(b - a, c - a);
  return 0.25 * (n * n);
}

inline double squared_distance(const Point3 & a, const Point3 & b)
{
  const Vector3 d = a - b;
  return d * d;
}

inline double to_double(double x) { return x; }
inline double sqrt(double x) { return std::sqrt(x); }

// orient3d: sign of the signed volume of (a,b,c,d), matching CGAL::orientation:
//   CGAL::orientation(a,b,c,d) == sign of det[(b-a),(c-a),(d-a)]
//                              == (d-a) . ((b-a) x (c-a)).
// > 0 (POSITIVE) when d is on the positive side of the oriented plane (a,b,c).
//
// TODO(Phase B): replace with Shewchuk's adaptive robust orient3d.
inline Sign orient3d(const Point3 & a, const Point3 & b, const Point3 & c,
                     const Point3 & d)
{
  const Vector3 ba = b - a;
  const Vector3 ca = c - a;
  const Vector3 da = d - a;
  const double det = da * cross_product(ba, ca);
  return det > 0 ? POSITIVE : (det < 0 ? NEGATIVE : ZERO);
}

// Named `orientation` to match the CGAL spelling used at call sites.
inline Sign orientation(const Point3 & a, const Point3 & b, const Point3 & c,
                        const Point3 & d)
{ return orient3d(a, b, c, d); }

}  // namespace nat
