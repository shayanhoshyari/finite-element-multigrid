#if !defined(GRID_ITERATORS_HXX_IS_INCLUDED)
#define GRID_ITERATORS_HXX_IS_INCLUDED

#include <cassert>
#include <vector>

#include "grid_connectivity_context.hxx"

// I feel like a librarian! So much book-keeping
// No support for combinatorial editing right no(delete, add, flip edge, etc.)
//
// TO CHECK: 1) Check the reset_boundary function of Iterator::Vertex_umbrella.
//              he and he.next() after reset must lie on the boundary
class Grid_iterators
{
public:
  class Half_edge;
  class Vertex;
  class Triangle;
  class Vertex_umbrella;

  Grid_iterators(const Grid_connectivity_context &);

  //
  const Grid_connectivity_context & connectivity() const { return *_connectivity_context; }

  // Give iterators
  int n_vertices() const;
  //
  int n_real_triangles() const;
  int n_dummy_triangles() const;
  int n_all_triangles() const;
  //
  int n_half_edges() const;
  int n_boundary_edges() const;
  int n_boundary_half_edges() const;
  int n_internal_half_edges() const;
  //
  Vertex vertex_iterator(const int) const;
  Triangle triangle_iterator(const int) const;
  Half_edge half_edge_iterator(const int, const int) const;
  Half_edge half_edge_iterator(const int) const;
  Half_edge boundary_half_edge_iterator(const int) const;
  Half_edge internal_half_edges_iterator(const int) const;
  Vertex_umbrella vertex_umbrella_iterator(const int) const;
  //
  // bool connecting_half_edge


  // Half-edge iterators
  class Vertex
  {
  public:
    Vertex() = default;
    Vertex(const Grid_connectivity_context &, const int id);
    Half_edge half_edge() const;
    Triangle triangle() const;
    int id() const;
    int tag() const;
    double x() const;
    double y() const;
    double z() const;
    const Grid_connectivity_context & manager() const { return *_manager; }
    bool is_boundary() const;
    bool is_boundary(int & tag0, int &tag1) const;

  private:
    // Using references(&) rather than pointers(*) sucks, since it deletes the default
    // constructor, copy constructor, etc. Also, it makes assignment operators a nightmare
    // if not impossible.
    Grid_connectivity_context const * _manager;
    int _vert_id;
  };

  class Half_edge
  {
  public:
    Half_edge() = default;
    Half_edge(const Grid_connectivity_context &, const int tri_id, const int offset);
    static Half_edge build_from_id(const Grid_connectivity_context &, const int id);
    static Half_edge build_from_internal_id(const Grid_connectivity_context &, const int id);
    static Half_edge build_from_boundary_id(const Grid_connectivity_context &, const int id);
    Half_edge twin() const;
    Half_edge next() const;
    Half_edge prev() const;
    Vertex origin() const;
    Triangle triangle() const;
    //
    int raw_half_edge_id() const;
    int specific_half_edge_id() const;
    int raw_edge_id() const;
    int specific_edge_id() const;
    bool is_boundary() const;
    //
    const Grid_connectivity_context & manager() const { return *_manager; }

  private:
    Grid_connectivity_context const * _manager;
    int _tri_id;
    int _offset;
  };

  class Triangle
  {
  public:
    Triangle() = default;
    Triangle(const Grid_connectivity_context &, const int);
    Half_edge half_edge() const;
    bool is_dummy() const;
    int tag() const;
    Vertex vertex_from_offset(const int) const;
    int has_vertex(const Vertex) const;
    Triangle neighbour_from_offset(const int) const;
    int has_neighbour(const Triangle) const;

    // Returns triangle id for internal triangles
    // Returns boundary edge id for boundary dummy triangle
    int specific_id() const;
    int raw_id() const;

    //
    const Grid_connectivity_context & manager() const { return *_manager; }

  private:
    Grid_connectivity_context const * _manager;
    int _tri_id;
  };

  class Vertex_umbrella
  {
  public:
    Vertex_umbrella() = default;
    // The Vertex_umbrella iterator iterates over the neighbours
    // of the vertex at the **HEAD** and not origin of current.
    Vertex_umbrella(Half_edge current);
    bool reset_boundary();
    bool advance();
    Half_edge half_edge() const { return _cur; }

  private:
    Half_edge _cur;
    Half_edge _end;
  };

private:
  Grid_connectivity_context const * _connectivity_context;
};

// Use inline functions in optimized mode
#if defined(NDEBUG)
#  include "grid_iterators.cxx"
#endif

#endif /* GRID_ITERATORS_HXX_IS_INCLUDED */
