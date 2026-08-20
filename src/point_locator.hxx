#if !defined(POINT_LOCATOR_HXX_IS_INCLUDED)
#define POINT_LOCATOR_HXX_IS_INCLUDED

#include <set>

#include "bridson/vec.h"
#include "fake_fprintf.hxx"
#include "grid_connectivity_context.hxx"

//
//
// CAVEAT 1: A vertex can be incident on two edges with different boundary tags.
//           For now this is ignored.
class Point_locator
{
public:
  typedef bridson::Vec3d Vec3d;
  typedef bridson::Vec2d Vec2d;

  enum
  {
    LOCATOR_SUCCESS = 0,
    LOCATOR_FAILURE,
    LOCATOR_FATAL_ERROR
  };

  Point_locator(const Grid_connectivity_context & conn);

  //
  // Look inside a single triangle
  //
  static int triangle_search(const Vec2d & x, const Vec2d & p0, const Vec2d & p1, const Vec2d & p2, Vec3d & xi);
  int triangle_search(const Vec2d & x, const int tri_id, Vec3d & xi) const;
  //
  static int edge_search(const Vec2d & x, const Vec2d & p0, const Vec2d & p1, double & xi, double & dist2);
  int boundary_edge_search(const Vec2d & x, const int bedge_id, double & xi, double & dist2) const;

  //
  // Locate a point inside the grid
  //
  struct Tri_query_in
  {
    Vec2d x;
    int tri_id;
    int max_bfs_layers;
    Tri_query_in() = default;
    Tri_query_in(const Tri_query_in &) = default;
    Tri_query_in(const Vec2d _x, const int _tri_id, const int _max_bfs);
  };
  struct Tri_query_out
  {
    Vec3d xi;
    int tri_id;
    Tri_query_out() = default;
    Tri_query_out(const Vec3d _xi, const int _tri_id);
  };
  //
  int directional_inside_search_helper(Tri_query_in, Tri_query_out &, int & n_attempts, std::set<int> & visited) const;
  int directional_inside_search(Tri_query_in, Tri_query_out &) const;
  int bfs_inside_search(Tri_query_in, Tri_query_out &) const;
  int bruteforce_inside_search(Tri_query_in, Tri_query_out &) const;

  //
  // Project a point on the boundary
  //
  struct Edge_query_in
  {
    Vec2d x;
    int tag;
    int bedge_id;
    int max_bfs_layers;
    Edge_query_in() = default;
    Edge_query_in(const Edge_query_in &) = default;
    Edge_query_in(const Vec2d _x, const int _tag, const int _boundary_edge, const int _max_bfs);
  };
  struct Edge_query_out
  {
    // I realized this inconsistency too late
    // The point location would be (1-xi)*v0+xi*v1
    double xi;
    int bedge_id;
    double dist2;
    Edge_query_out() = default;
    Edge_query_out(const double _xi, const int _edge, const double _dist2);
  };
  //
  int bfs_boundary_search(Edge_query_in, Edge_query_out &) const;
  int bruteforce_boundary_search(Edge_query_in, Edge_query_out &) const;

  //
  // Project one mesh on another one
  //
  struct Grid_query_out
  {
    // Triangle id of the owner. Will be a ghost triangle for boundary
    // vertices.
    std::vector<int> owner_id;
    // Barycentric or normal coordinates. With block size of 2.
    // Use 1- xi[0] - xi[1] to find the third one.
    std::vector<double> xi;
  };
  //
  // FULL OF BUGS, all the umbrella iterators are wrong
  //
  void bfs_mesh_boundary_propagative_search(const Grid_connectivity_context & inconn,
                                            const int seed_vertex_id,
                                            Point_locator::Grid_query_out & oo) const;
  void bfs_mesh_internal_propagative_search(const Grid_connectivity_context & inconn,
                                            const int seed_vertex_id,
                                            const int seed_mytri_guess,
                                            Point_locator::Grid_query_out & oo) const;
  int bfs_mesh_search(const Grid_connectivity_context & in, Grid_query_out & oo) const;
  //
  int bruteforce_mesh_search(const Grid_connectivity_context & in, Grid_query_out & oo) const;

  //
  const Grid_connectivity_context & connectivity() const { return *_connectivity; }
  static void turn_on_debugging();
  static void turn_off_debugging();
  static constexpr double tol() { return 1e-10; }

private:
  bool is_valid_btag(const int i) const
  {
    return std::find(connectivity().boundary_tags().cbegin(), connectivity().boundary_tags().cend(), i) !=
           connectivity().boundary_tags().cend();
  };
  Grid_connectivity_context const * _connectivity;
  // Naivest possible way of debugging
  static Fprintf_function _printf_functor;
  // Do something about this later
  // No need to worry about a sqaure mesh however.
  constexpr static int MAX_BOUNDARY_BFS_LAYERS = 100;
};


#endif /* POINT_LOCATOR_HXX_IS_INCLUDED */