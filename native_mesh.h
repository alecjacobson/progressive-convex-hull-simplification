#pragma once
// -----------------------------------------------------------------------------
// Native half-edge surface mesh (Backend B) — std-only, index-based.
//
// Design goals (informed by what was painful in CGAL):
//   * The operation that "copies a candidate retriangulation into the main
//     mesh" is a FIRST-CLASS, easy primitive: retriangulate_star(v, tris)
//     erases v's triangle fan and fills the resulting hole directly from a list
//     of triangles over the ring vertices. No ear-clipping path to record and
//     replay.
//   * Opposite half-edges are paired: opposite(i) == (i ^ 1). Edges are just
//     even/odd half-edge pairs, so edge iteration and twin lookup are trivial.
//   * Handles are lightweight {mesh, index} values that also act as forward
//     iterators (operator++ skips deleted slots), mirroring how CGAL handles are
//     iterators. operator-> returns the handle itself, so `h->next()` works.
//
// Border half-edges are represented explicitly (facet index = -1), so open
// meshes (candidate one-rings, holes mid-operation) are first class.
//
// Deletion marks slots and appends new ones; call garbage_collect() to compact.
// -----------------------------------------------------------------------------
#include "native_geom.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace nat
{

class Mesh;
class EdgeIter;

// Shared handle: a {mesh, index} pair that also serves as a forward iterator.
// Kind selects which array it indexes (for iteration bounds / skip-deleted).
enum class Elem { V, H, F };

template <Elem K>
class Handle
{
public:
  Handle() = default;
  Handle(Mesh * m, int i) : m_(m), i_(i) {}

  // Acts as its own element proxy: `h->method()` == `h.method()`.
  const Handle * operator->() const { return this; }
  const Handle & operator*() const { return *this; }

  bool operator==(const Handle & o) const { return m_ == o.m_ && i_ == o.i_; }
  bool operator!=(const Handle & o) const { return !(*this == o); }
  bool operator==(std::nullptr_t) const { return i_ < 0; }
  bool operator!=(std::nullptr_t) const { return i_ >= 0; }
  // Total order over the index (the flip loop relies on comparing half-edges).
  bool operator<(const Handle & o) const { return i_ < o.i_; }
  bool operator>(const Handle & o) const { return i_ > o.i_; }

  Handle & operator++();                 // defined after Mesh
  int index() const { return i_; }
  Mesh * mesh() const { return m_; }

  // --- element accessors (only the valid ones for K are defined below) ---
  // Vertex:
  Handle<Elem::H> halfedge() const;      // K == V
  int & id() const;                      // all
  Point3 & point() const;                // K == V
  int degree() const;                    // K == V
  // Halfedge:
  Handle<Elem::V> vertex() const;        // K == H (target vertex)
  Handle<Elem::H> next() const;          // K == H
  Handle<Elem::H> prev() const;          // K == H
  Handle<Elem::H> opposite() const;      // K == H
  Handle<Elem::F> facet() const;         // K == H
  bool is_border() const;                // K == H
  bool is_border_edge() const;           // K == H
  // Facet:
  Handle<Elem::H> facet_halfedge() const;// K == F  (spelled halfedge() too)

private:
  Mesh * m_ = nullptr;
  int i_ = -1;
};

class Mesh
{
public:
  using Point_3 = Point3;
  struct Traits { using FT = double; using Point_3 = Point3;
                  using Vector_3 = Vector3; using Oriented_side = Sign; };

  using Vertex_handle          = Handle<Elem::V>;
  using Vertex_const_handle    = Handle<Elem::V>;
  using Vertex_iterator        = Handle<Elem::V>;
  using Vertex_const_iterator  = Handle<Elem::V>;
  using Halfedge_handle        = Handle<Elem::H>;
  using Halfedge_const_handle  = Handle<Elem::H>;
  using Halfedge_iterator      = Handle<Elem::H>;
  using Halfedge_const_iterator= Handle<Elem::H>;
  using Facet_handle           = Handle<Elem::F>;
  using Facet_const_handle     = Handle<Elem::F>;
  using Facet_iterator         = Handle<Elem::F>;
  using Facet_const_iterator   = Handle<Elem::F>;

  struct HE { int next = -1, prev = -1, v = -1, f = -1; int id = 0; bool dead = false; };
  struct VR { int he = -1; int id = 0; Point3 p; bool dead = false; };
  struct FR { int he = -1; int id = 0; bool dead = false; };

  // ---- storage (public so Handle can read; treat as private) ----
  std::vector<HE> he_;
  std::vector<VR> vr_;
  std::vector<FR> fr_;
  int nv_ = 0, nh_ = 0, nf_ = 0;   // live counts

  static constexpr int NULLID = -1;

  Mesh() = default;

  // --- sizes / predicates ---
  std::size_t size_of_vertices()  const { return nv_; }
  std::size_t size_of_halfedges() const { return nh_; }
  std::size_t size_of_facets()    const { return nf_; }

  bool is_closed() const
  {
    for(int i = 0; i < (int)he_.size(); ++i)
      if(!he_[i].dead && he_[i].f < 0) return false;
    return true;
  }
  bool is_pure_triangle() const
  {
    for(int f = 0; f < (int)fr_.size(); ++f)
    {
      if(fr_[f].dead) continue;
      int h = fr_[f].he, c = 0;
      do { h = he_[h].next; ++c; } while(h != fr_[f].he && c < 1000);
      if(c != 3) return false;
    }
    return true;
  }
  // Euler characteristic V - E + F (E = live halfedges / 2).
  int euler() const { return nv_ - nh_ / 2 + nf_; }

  void reserve(std::size_t v, std::size_t h, std::size_t f)
  { vr_.reserve(v); he_.reserve(h); fr_.reserve(f); }

  // --- iteration (skip deleted) ---
  template <Elem K> Handle<K> begin_() const
  {
    Mesh * m = const_cast<Mesh *>(this);
    Handle<K> h(m, -1);
    return ++h;  // advance to first live (operator++ handles -1 -> 0 scan)
  }
  Vertex_iterator  vertices_begin()  const { return begin_<Elem::V>(); }
  Vertex_iterator  vertices_end()    const { return {const_cast<Mesh*>(this), (int)vr_.size()}; }
  Halfedge_iterator halfedges_begin() const { return begin_<Elem::H>(); }
  Halfedge_iterator halfedges_end()   const { return {const_cast<Mesh*>(this), (int)he_.size()}; }
  Facet_iterator   facets_begin()    const { return begin_<Elem::F>(); }
  Facet_iterator   facets_end()      const { return {const_cast<Mesh*>(this), (int)fr_.size()}; }
  // Edge iteration: one representative half-edge per edge. Steps by 2 (opposite
  // pairs are adjacent), so each edge is visited exactly once.
  EdgeIter edges_begin() const;
  EdgeIter edges_end()   const;

  // --- construction ------------------------------------------------------
  // Build from a triangle soup. `tris` are CCW index triples into `pts`.
  // Open edges (used by only one triangle) get border half-edges.
  void build(const std::vector<Point3> & pts,
             const std::vector<std::array<int, 3>> & tris);

  // Fill hole / retriangulate: erase vertex v's triangle fan and fill the
  // resulting ring polygon with `tris` (triples of indices into the ring order
  // returned by ring_vertices(v)). THE easy splice.
  void retriangulate_star(Vertex_handle v,
                          const std::vector<std::array<int, 3>> & tris);

  // Ordered ring (link) of v: neighbor vertex handles around v.
  std::vector<Vertex_handle> ring_vertices(Vertex_handle v) const;

  // Flip the interior edge of half-edge h (both incident facets are triangles).
  void flip_edge(Halfedge_handle h);

  // Compact storage, dropping deleted slots. Invalidates handles/indices.
  void garbage_collect();

  // --- low-level slot allocation (used by build/retriangulate) ---
  int new_vertex(const Point3 & p) { vr_.push_back({-1, 0, p, false}); ++nv_; return (int)vr_.size()-1; }
  int new_facet() { fr_.push_back({-1, 0, false}); ++nf_; return (int)fr_.size()-1; }
  // Allocate an opposite half-edge pair; returns the even index e (odd is e^1).
  int new_edge() { int e = (int)he_.size(); he_.push_back({}); he_.push_back({}); nh_ += 2; return e; }
  void kill_halfedge(int i) { if(!he_[i].dead){ he_[i].dead = true; --nh_; } }
  void kill_vertex(int i)   { if(!vr_[i].dead){ vr_[i].dead = true; --nv_; } }
  void kill_facet(int i)    { if(!fr_[i].dead){ fr_[i].dead = true; --nf_; } }
};

// ---- Handle method definitions (need full Mesh) ----

template <Elem K>
inline Handle<K> & Handle<K>::operator++()
{
  const int n = (K == Elem::V) ? (int)m_->vr_.size()
              : (K == Elem::H) ? (int)m_->he_.size()
              :                  (int)m_->fr_.size();
  ++i_;
  while(i_ < n)
  {
    const bool dead = (K == Elem::V) ? m_->vr_[i_].dead
                    : (K == Elem::H) ? m_->he_[i_].dead
                    :                  m_->fr_[i_].dead;
    if(!dead) break;
    ++i_;
  }
  return *this;
}

template <> inline Handle<Elem::H> Handle<Elem::V>::halfedge() const
{ return {m_, m_->vr_[i_].he}; }
template <> inline Point3 & Handle<Elem::V>::point() const
{ return m_->vr_[i_].p; }
template <> inline int & Handle<Elem::V>::id() const
{ return m_->vr_[i_].id; }
template <> inline int & Handle<Elem::H>::id() const
{ return m_->he_[i_].id; }
template <> inline int & Handle<Elem::F>::id() const
{ return m_->fr_[i_].id; }

template <> inline Handle<Elem::V> Handle<Elem::H>::vertex() const
{ return {m_, m_->he_[i_].v}; }
template <> inline Handle<Elem::H> Handle<Elem::H>::next() const
{ return {m_, m_->he_[i_].next}; }
template <> inline Handle<Elem::H> Handle<Elem::H>::prev() const
{ return {m_, m_->he_[i_].prev}; }
template <> inline Handle<Elem::H> Handle<Elem::H>::opposite() const
{ return {m_, i_ ^ 1}; }
template <> inline Handle<Elem::F> Handle<Elem::H>::facet() const
{ return {m_, m_->he_[i_].f}; }
template <> inline bool Handle<Elem::H>::is_border() const
{ return m_->he_[i_].f < 0; }
template <> inline bool Handle<Elem::H>::is_border_edge() const
{ return m_->he_[i_].f < 0 || m_->he_[i_ ^ 1].f < 0; }
template <> inline Handle<Elem::H> Handle<Elem::F>::halfedge() const
{ return {m_, m_->fr_[i_].he}; }

template <> inline int Handle<Elem::V>::degree() const
{
  // v->halfedge() is incoming (target v); rotate via next()->opposite().
  const int start = m_->vr_[i_].he;
  if(start < 0) return 0;
  int d = 0, h = start;
  do { ++d; h = m_->he_[h].next ^ 1; } while(h != start && d < 1000000);
  return d;
}

// One representative half-edge per edge; operator++ steps by 2 (skipping dead
// pairs). Dereferences to a Halfedge handle so callers use e->next() etc.
class EdgeIter
{
public:
  EdgeIter(Mesh * m, int i) : m_(m), i_(i) {}
  Handle<Elem::H> operator*() const { return {m_, i_}; }
  Handle<Elem::H> operator->() const { return {m_, i_}; }
  bool operator==(const EdgeIter & o) const { return i_ == o.i_; }
  bool operator!=(const EdgeIter & o) const { return i_ != o.i_; }
  EdgeIter & operator++()
  {
    const int n = (int)m_->he_.size();
    i_ += 2;
    while(i_ < n && m_->he_[i_].dead) i_ += 2;
    return *this;
  }
private:
  Mesh * m_;
  int i_;
};

inline EdgeIter Mesh::edges_begin() const
{
  Mesh * m = const_cast<Mesh *>(this);
  int i = 0;
  while(i < (int)he_.size() && he_[i].dead) i += 2;
  return {m, i};
}
inline EdgeIter Mesh::edges_end() const
{
  return {const_cast<Mesh *>(this), (int)he_.size()};
}

// ---------------------------------------------------------------------------
// build: triangle soup -> half-edge mesh (border half-edges for open edges)
// ---------------------------------------------------------------------------
inline void Mesh::build(const std::vector<Point3> & pts,
                        const std::vector<std::array<int, 3>> & tris)
{
  he_.clear(); vr_.clear(); fr_.clear(); nv_ = nh_ = nf_ = 0;
  for(const auto & p : pts) new_vertex(p);

  // Undirected edge {min,max} -> even half-edge index (dir min->max is e).
  std::map<std::pair<int, int>, int> edge;
  auto half_for = [&](int u, int w) -> int
  {
    const std::pair<int, int> key(std::min(u, w), std::max(u, w));
    auto it = edge.find(key);
    int e;
    if(it == edge.end())
    {
      e = new_edge();
      he_[e].v     = key.second;   // dir key.first -> key.second
      he_[e ^ 1].v = key.first;
      edge.emplace(key, e);
    }
    else e = it->second;
    return (u < w) ? e : (e ^ 1);   // half-edge directed u->w
  };

  for(const auto & t : tris)
  {
    const int fid = new_facet();
    int h[3];
    for(int k = 0; k < 3; ++k) h[k] = half_for(t[k], t[(k + 1) % 3]);
    for(int k = 0; k < 3; ++k)
    {
      he_[h[k]].f    = fid;
      he_[h[k]].next = h[(k + 1) % 3];
      he_[h[k]].prev = h[(k + 2) % 3];
      vr_[t[(k + 1) % 3]].he = h[k];        // incoming interior half-edge
    }
    fr_[fid].he = h[0];
  }

  // Link border half-edges (facet still -1) into boundary loops.
  std::map<int, int> out_of;   // source vertex -> border half-edge
  for(int b = 0; b < (int)he_.size(); ++b)
    if(!he_[b].dead && he_[b].f < 0)
      out_of[he_[b ^ 1].v] = b;
  for(int b = 0; b < (int)he_.size(); ++b)
    if(!he_[b].dead && he_[b].f < 0)
    {
      const int nb = out_of[he_[b].v];       // next border starts at b's target
      he_[b].next = nb;
      he_[nb].prev = b;
    }
}

// ---------------------------------------------------------------------------
// boundary_of_star: ordered ring boundary half-edges B[] and ring vertices p[]
// (global indices). B[j] is directed p[(j-1+k)%k] -> p[j] and its opposite is
// an unaffected outside face; used by ring_vertices and retriangulate_star.
// ---------------------------------------------------------------------------
inline void boundary_of_star(Mesh & M, Mesh::Vertex_handle v,
                             std::vector<int> & B, std::vector<int> & p)
{
  B.clear(); p.clear();
  const int start = M.vr_[v.index()].he;   // incoming, target v
  // Star-side boundary half-edges: third edge of each incident triangle.
  std::vector<int> bh;
  std::map<int, int> src_to_bh;            // source vertex -> boundary half-edge
  int h = start;
  do {
    const int b = M.he_[M.he_[h].next].next;   // h->next->next
    bh.push_back(b);
    src_to_bh[M.he_[b ^ 1].v] = b;             // source(b) -> b
    h = M.he_[h].next ^ 1;                      // rotate around v
  } while(h != start);
  // Walk the boundary cycle in half-edge direction.
  int cur = bh[0];
  for(std::size_t i = 0; i < bh.size(); ++i)
  {
    B.push_back(cur);
    p.push_back(M.he_[cur].v);              // target
    cur = src_to_bh[M.he_[cur].v];          // next boundary edge starts here
  }
}

inline std::vector<Mesh::Vertex_handle> Mesh::ring_vertices(Vertex_handle v) const
{
  Mesh * m = const_cast<Mesh *>(this);
  std::vector<int> B, p;
  boundary_of_star(*m, v, B, p);
  std::vector<Vertex_handle> out;
  for(int gi : p) out.push_back({m, gi});
  return out;
}

// ---------------------------------------------------------------------------
// flip_edge: flip the interior edge of h (both incident facets are triangles)
// ---------------------------------------------------------------------------
inline void Mesh::flip_edge(Halfedge_handle hh)
{
  const int h = hh.index(), o = h ^ 1;
  assert(he_[h].f >= 0 && he_[o].f >= 0);
  const int f1 = he_[h].f, f2 = he_[o].f;
  const int n1 = he_[h].next, n2 = he_[n1].next;   // f1 = h,n1,n2
  const int m1 = he_[o].next, m2 = he_[m1].next;   // f2 = o,m1,m2
  const int v0 = he_[o].v, v1 = he_[h].v;
  const int v2 = he_[n1].v, v3 = he_[m1].v;

  // New diagonal v3<->v2: h becomes v3->v2, o becomes v2->v3.
  he_[h].v = v2; he_[o].v = v3;
  // f1 = (m1: v0->v3, h: v3->v2, n2: v2->v0)
  he_[m1].f = f1; he_[h].f = f1; he_[n2].f = f1;
  he_[m1].next = h; he_[h].next = n2; he_[n2].next = m1;
  he_[h].prev = m1; he_[n2].prev = h; he_[m1].prev = n2;
  // f2 = (n1: v1->v2, o: v2->v3, m2: v3->v1)
  he_[n1].f = f2; he_[o].f = f2; he_[m2].f = f2;
  he_[n1].next = o; he_[o].next = m2; he_[m2].next = n1;
  he_[o].prev = n1; he_[m2].prev = o; he_[n1].prev = m2;

  fr_[f1].he = h; fr_[f2].he = o;
  vr_[v0].he = n2; vr_[v1].he = m2; vr_[v2].he = n1; vr_[v3].he = m1;
}

// ---------------------------------------------------------------------------
// retriangulate_star: erase v's fan and fill the hole from `tris` (triples over
// ring_vertices(v) order). THE easy splice — no ear-clip path replay.
// ---------------------------------------------------------------------------
inline void Mesh::retriangulate_star(Vertex_handle v,
                                     const std::vector<std::array<int, 3>> & tris)
{
  std::vector<int> B, p;
  boundary_of_star(*this, v, B, p);
  const int k = (int)p.size();

  // Delete v, its spoke edges, and the star facets. Keep the B[] boundary
  // half-edges (their opposites are unaffected outside faces).
  {
    const int start = vr_[v.index()].he;
    int h = start;
    do {
      kill_facet(he_[h].f);          // star facet incident to h
      const int nh = he_[h].next ^ 1;
      kill_halfedge(h); kill_halfedge(h ^ 1);   // spoke edge (v, source(h))
      h = nh;
    } while(h != start);
    kill_vertex(v.index());
  }

  // Local ring index -> global vertex.
  auto gv = [&](int a) { return p[a]; };

  // Half-edge for ring edge a->b: boundary edge reuses B, diagonals are new.
  std::map<std::pair<int, int>, int> diag;
  auto he_for = [&](int a, int b) -> int
  {
    if(b == (a + 1) % k) return B[b];        // boundary p_a -> p_b (== B[b])
    auto key = std::make_pair(a, b);
    auto it = diag.find(key);
    if(it != diag.end()) return it->second;
    const int e = new_edge();
    he_[e].v = gv(b); he_[e ^ 1].v = gv(a);
    diag[{a, b}] = e; diag[{b, a}] = e ^ 1;
    return e;
  };

  for(const auto & t : tris)
  {
    const int fid = new_facet();
    int hh[3];
    for(int j = 0; j < 3; ++j) hh[j] = he_for(t[j], t[(j + 1) % 3]);
    for(int j = 0; j < 3; ++j)
    {
      he_[hh[j]].f = fid;
      he_[hh[j]].next = hh[(j + 1) % 3];
      he_[hh[j]].prev = hh[(j + 2) % 3];
    }
    fr_[fid].he = hh[0];
  }
  for(int j = 0; j < k; ++j) vr_[gv(j)].he = B[j];   // incoming interior he
}

// ---------------------------------------------------------------------------
// garbage_collect: compact arrays, dropping deleted slots. Reassigns indices.
// ---------------------------------------------------------------------------
inline void Mesh::garbage_collect()
{
  std::vector<int> vmap(vr_.size(), -1), fmap(fr_.size(), -1), hmap(he_.size(), -1);
  std::vector<VR> nvr; std::vector<FR> nfr; std::vector<HE> nhe;
  for(int i = 0; i < (int)vr_.size(); ++i)
    if(!vr_[i].dead) { vmap[i] = (int)nvr.size(); nvr.push_back(vr_[i]); }
  for(int i = 0; i < (int)fr_.size(); ++i)
    if(!fr_[i].dead) { fmap[i] = (int)nfr.size(); nfr.push_back(fr_[i]); }
  // Keep edge pairs together so opposite(i)==i^1 survives compaction.
  for(int e = 0; e < (int)he_.size(); e += 2)
    if(!he_[e].dead || !he_[e + 1].dead)
    {
      hmap[e] = (int)nhe.size(); hmap[e + 1] = (int)nhe.size() + 1;
      nhe.push_back(he_[e]); nhe.push_back(he_[e + 1]);
    }
  auto rh = [&](int i) { return i < 0 ? -1 : hmap[i]; };
  for(auto & h : nhe) { h.next = rh(h.next); h.prev = rh(h.prev);
                        h.v = h.v < 0 ? -1 : vmap[h.v];
                        h.f = h.f < 0 ? -1 : fmap[h.f]; }
  for(auto & v : nvr) v.he = rh(v.he);
  for(auto & f : nfr) f.he = rh(f.he);
  vr_ = std::move(nvr); fr_ = std::move(nfr); he_ = std::move(nhe);
  nv_ = (int)vr_.size(); nf_ = (int)fr_.size();
  nh_ = 0; for(auto & h : he_) if(!h.dead) ++nh_;
}

}  // namespace nat
