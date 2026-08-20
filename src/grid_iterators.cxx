#if defined(NDEBUG)
#define GITER_INLINE inline
#else
#include "grid_iterators.hxx"
#define GITER_INLINE
#endif

namespace giter_static
{
constexpr static int plus_1[] = {1, 2, 0};
constexpr static int minus_1[] = {2, 0, 1};
}

// ===================================================
//            GRID
// ===================================================
GITER_INLINE
Grid_iterators::Grid_iterators(const Grid_connectivity_context & d)
: _connectivity_context(&d)
{
}


// ======================= Iterators
#define DUPLICATE_ME(name, rtype) \
  rtype GITER_INLINE Grid_iterators::name() const { return connectivity().name(); }

DUPLICATE_ME(n_vertices, int)
DUPLICATE_ME(n_real_triangles, int)
DUPLICATE_ME(n_dummy_triangles, int)
DUPLICATE_ME(n_all_triangles, int)

#undef DUPLICATE_ME


GITER_INLINE int
Grid_iterators::n_half_edges() const
{
  return n_boundary_half_edges() + n_internal_half_edges();
}

GITER_INLINE int
Grid_iterators::n_boundary_half_edges() const
{
  return n_dummy_triangles();
}

GITER_INLINE int
Grid_iterators::n_internal_half_edges() const
{
  return n_real_triangles() * 3;
}

GITER_INLINE int
Grid_iterators::n_boundary_edges() const
{
  return n_boundary_half_edges();
}


GITER_INLINE Grid_iterators::Vertex
Grid_iterators::vertex_iterator(const int vid) const
{
  return Vertex(this->connectivity(), vid);
}

GITER_INLINE Grid_iterators::Triangle
Grid_iterators::triangle_iterator(const int tid) const
{
  return Triangle(this->connectivity(), tid);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::half_edge_iterator(const int tri_id, const int offset) const
{
  return Half_edge(this->connectivity(), tri_id, offset);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::half_edge_iterator(const int hid) const
{
  return Half_edge::build_from_id(this->connectivity(), hid);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::boundary_half_edge_iterator(const int hid) const
{
  return Half_edge::build_from_boundary_id(this->connectivity(), hid);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::internal_half_edges_iterator(const int hid) const
{
  return Half_edge::build_from_internal_id(this->connectivity(), hid);
}

GITER_INLINE Grid_iterators::Vertex_umbrella
Grid_iterators::vertex_umbrella_iterator(const int vid) const
{
  return Vertex_umbrella(vertex_iterator(vid).half_edge().twin());
}


// ===================================================
//         VERTEX ITERATOR (book-keeping is here)
// ===================================================

GITER_INLINE
Grid_iterators::Vertex::Vertex(const Grid_connectivity_context & g, const int id)
: _manager(&g)
, _vert_id(id)
{
  assert(_vert_id < manager().n_vertices());
}

GITER_INLINE Grid_iterators::Triangle
Grid_iterators::Vertex::triangle() const
{
  // works on boundary?
  const Grid_connectivity_context::Vertex vertex_data = manager().vertices_data()[_vert_id];
  return Grid_iterators(manager()).triangle_iterator(vertex_data.tri_id);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Vertex::half_edge() const
{
  Triangle tri = triangle();
  const int my_offset = tri.has_vertex(*this);
  assert(my_offset >= 0);
  assert(my_offset <= 2);
  return Grid_iterators(manager()).half_edge_iterator(tri.raw_id(), giter_static::minus_1[my_offset]);
}

GITER_INLINE int
Grid_iterators::Vertex::id() const
{
  return _vert_id;
}

GITER_INLINE int
Grid_iterators::Vertex::tag() const
{
  return manager().vertices_data()[_vert_id].tag;
}

GITER_INLINE double
Grid_iterators::Vertex::x() const
{
  return manager().vertices_data()[_vert_id].xy[0];
}

GITER_INLINE double
Grid_iterators::Vertex::y() const
{
  return manager().vertices_data()[_vert_id].xy[1];
}

GITER_INLINE double
Grid_iterators::Vertex::z() const
{
  return manager().vertices_data()[_vert_id].xy[2];
}

GITER_INLINE bool
Grid_iterators::Vertex::is_boundary() const
{
  return Grid_iterators(manager()).vertex_umbrella_iterator(this->id()).reset_boundary();
}

GITER_INLINE bool
Grid_iterators::Vertex::is_boundary(int & tag0, int & tag1) const
{
  Grid_iterators::Vertex_umbrella umbrella_iter = Grid_iterators(manager()).vertex_umbrella_iterator(this->id());
  const bool answer = umbrella_iter.reset_boundary();
  if(answer)
    {
      tag0 = umbrella_iter.half_edge().triangle().tag();
      tag1 = umbrella_iter.half_edge().next().triangle().tag();
      assert(umbrella_iter.half_edge().next().triangle().is_dummy());
    }
  return answer;
}

// ===================================================
//         TRIANGLE ITERATOR (book-keeping is here)
// ===================================================


GITER_INLINE
Grid_iterators::Triangle::Triangle(const Grid_connectivity_context & g, const int tid)
: _manager(&g)
, _tri_id(tid)
{
  assert(_tri_id < manager().n_all_triangles());
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Triangle::half_edge() const
{
  // 0 will work for both real and dummy.
  // 1 and 2 will only work for internal.
  const int vertex_offset = 0;
  return Grid_iterators(manager()).half_edge_iterator(_tri_id, vertex_offset);
}

GITER_INLINE bool
Grid_iterators::Triangle::is_dummy() const
{
  return _tri_id >= manager().n_real_triangles();
}

GITER_INLINE int
Grid_iterators::Triangle::tag() const
{
  return manager().triangles_data()[_tri_id].tag;
}

GITER_INLINE Grid_iterators::Vertex
Grid_iterators::Triangle::vertex_from_offset(const int offset) const
{
  assert(offset >= 0);
  assert((!is_dummy()) || offset >= 1);
  assert(offset <= 2);
  const int vertex_id = manager().triangles_data()[_tri_id].vert_ids[offset];
  return Grid_iterators(manager()).vertex_iterator(vertex_id);
}

GITER_INLINE int
Grid_iterators::Triangle::has_vertex(const Vertex vin) const
{
  if(vertex_from_offset(1).id() == vin.id())
    return 1;
  else if(vertex_from_offset(2).id() == vin.id())
    return 2;
  else if((!is_dummy()) && (vertex_from_offset(0).id() == vin.id()))
    return 0;
  else
    return -1;
}

GITER_INLINE Grid_iterators::Triangle
Grid_iterators::Triangle::neighbour_from_offset(const int offset) const
{
  assert(offset >= 0);
  assert(offset <= 2);
  const int neighbour_id = manager().triangles_data()[_tri_id].tri_ids[offset];
  return Grid_iterators(manager()).triangle_iterator(neighbour_id);
}

GITER_INLINE int
Grid_iterators::Triangle::has_neighbour(const Triangle tin) const
{
  if(neighbour_from_offset(0).raw_id() == tin.raw_id())
    return 0;
  else if((!is_dummy()) && (neighbour_from_offset(1).raw_id() == tin.raw_id()))
    return 1;
  else if((!is_dummy()) && (neighbour_from_offset(2).raw_id() == tin.raw_id()))
    return 2;
  else
    return -1;
}

GITER_INLINE int
Grid_iterators::Triangle::specific_id() const
{
  if(is_dummy())
    return raw_id() - manager().n_real_triangles();
  else
    return raw_id();
}

GITER_INLINE int
Grid_iterators::Triangle::raw_id() const
{
  return _tri_id;
}

// ===================================================
//         HALF-EDGE ITERATOR (book-keeping is here)
// ===================================================


GITER_INLINE
Grid_iterators::Half_edge::Half_edge(const Grid_connectivity_context & g, const int tri_id, const int offset)
: _manager(&g)
, _tri_id(tri_id)
, _offset(offset)
{
  assert((!triangle().is_dummy()) || (_offset == 0));
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Half_edge::build_from_id(const Grid_connectivity_context & c, const int id)
{
  Grid_iterators g(c);
  assert(id < g.n_half_edges());

  if(id < g.n_internal_half_edges())
    {
      const int tri_id = id / 3;
      const int offset = id % 3;
      return Half_edge(c, tri_id, offset);
    }
  else
    {
      const int tri_id = g.n_real_triangles() + (id - g.n_internal_half_edges());
      const int offset = 0;
      return Half_edge(c, tri_id, offset);
    }
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Half_edge::build_from_internal_id(const Grid_connectivity_context & g, const int id)
{
  assert(id < Grid_iterators(g).n_internal_half_edges());
  return build_from_id(g, id);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Half_edge::build_from_boundary_id(const Grid_connectivity_context & g, const int id)
{
  assert(id < Grid_iterators(g).n_boundary_half_edges());
  return build_from_id(g, id + Grid_iterators(g).n_internal_half_edges());
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Half_edge::twin() const
{
  const Triangle tri = triangle();
  const Triangle tri_n = tri.neighbour_from_offset(_offset);
  const int my_origins_offset_on_neigh = tri_n.has_vertex(origin());
  assert(my_origins_offset_on_neigh >= 0);
  return Grid_iterators(manager()).half_edge_iterator(tri_n.raw_id(), giter_static::plus_1[my_origins_offset_on_neigh]);
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Half_edge::next() const
{
  if(is_boundary())
    {
      const Triangle tri = triangle();
      //
      assert(tri.is_dummy());
      assert(_offset == 0);
      const Triangle tri_n = tri.neighbour_from_offset(giter_static::plus_1[_offset]);
      assert(tri_n.is_dummy());
      //
      const Vertex v_common = tri.vertex_from_offset(giter_static::minus_1[_offset]);
      const int v_common_on_neighbour_offset = tri_n.has_vertex(v_common);
      assert(v_common_on_neighbour_offset == 1);
      //
      return Grid_iterators(manager()).half_edge_iterator(tri_n.raw_id(),
                                                          giter_static::minus_1[v_common_on_neighbour_offset]);
    }
  else
    {
      return Grid_iterators(manager()).half_edge_iterator(_tri_id, giter_static::plus_1[_offset]);
    }
}

GITER_INLINE Grid_iterators::Half_edge
Grid_iterators::Half_edge::prev() const
{
  if(is_boundary())
    {
      const Triangle tri = triangle();
      //
      assert(tri.is_dummy());
      assert(_offset == 0);
      const Triangle tri_n = tri.neighbour_from_offset(giter_static::minus_1[_offset]);
      assert(tri_n.is_dummy());
      //
      const Vertex v_common = tri.vertex_from_offset(giter_static::plus_1[_offset]);
      const int v_common_on_neighbour_offset = tri_n.has_vertex(v_common);
      assert(v_common_on_neighbour_offset == 2);
      //
      return Grid_iterators(manager()).half_edge_iterator(tri_n.raw_id(),
                                                          giter_static::plus_1[v_common_on_neighbour_offset]);
    }
  else
    {
      return Grid_iterators(manager()).half_edge_iterator(_tri_id, giter_static::minus_1[_offset]);
    }
}

GITER_INLINE Grid_iterators::Vertex
Grid_iterators::Half_edge::origin() const
{
  return triangle().vertex_from_offset(giter_static::plus_1[_offset]);
}

GITER_INLINE Grid_iterators::Triangle
Grid_iterators::Half_edge::triangle() const
{
  return Grid_iterators(manager()).triangle_iterator(_tri_id);
}

GITER_INLINE int
Grid_iterators::Half_edge::specific_half_edge_id() const
{
  if(is_boundary())
    {
      Grid_iterators c(manager());
      return _tri_id - c.n_real_triangles();
    }
  else
    {
      return 3 * _tri_id + _offset;
    }
}

GITER_INLINE int
Grid_iterators::Half_edge::raw_half_edge_id() const
{
  if(is_boundary())
    {
      Grid_iterators c(manager());
      return _tri_id - c.n_real_triangles() + c.n_internal_half_edges();
    }
  else
    {
      return 3 * _tri_id + _offset;
    }
}


GITER_INLINE int
Grid_iterators::Half_edge::raw_edge_id() const
{
  // if(is_boundary())
  //  {
  //    return Grid_iterators(manager()).n_boundary_half_edges() / 2 + specific_half_edge_id();
  //  }
  // else
  {
    assert(0 && "not implemented; no edge book-keeping yet");
    throw;
  }
  return -1;
}

GITER_INLINE int
Grid_iterators::Half_edge::specific_edge_id() const
{
  // if(is_boundary())
  //  {
  //    return specific_half_edge_id();
  //  }
  // else
  {
    assert(0 && "not implemented; no edge book-keeping yet");
    throw;
  }
  return -1;
}

GITER_INLINE bool
Grid_iterators::Half_edge::is_boundary() const
{
  return triangle().is_dummy();
}

// ===================================================
//         VERTEX UMBRELLA ITERATOR (book-keeping is here)
// ===================================================
GITER_INLINE
Grid_iterators::Vertex_umbrella::Vertex_umbrella(Half_edge current)
: _cur(current)
, _end(current)
{
}

GITER_INLINE bool
Grid_iterators::Vertex_umbrella::reset_boundary()
{
  Half_edge snapshot = _cur;
  do
    {
      // start from the boundary edge located to the right of the vertex.
      if((_cur.is_boundary()) && (_cur.next().is_boundary()))
        {
          _end = _cur;
          return true;
        }
      _cur = _cur.next().twin();
    }
  while(_cur.raw_half_edge_id() != snapshot.raw_half_edge_id());
  return false;
}


GITER_INLINE bool
Grid_iterators::Vertex_umbrella::advance()
{
  _cur = _cur.next().twin();
  return _cur.raw_half_edge_id() != _end.raw_half_edge_id();
}


#undef GITER_INLINE