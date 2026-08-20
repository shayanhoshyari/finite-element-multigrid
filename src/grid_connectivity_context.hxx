#if !defined(GRID_CONNECTIVITY_CONTEXT)
#define GRID_CONNECTIVITY_CONTEXT

#include <cassert>
#include <vector>

#include "removable_object.hxx"

// TODOS: 1) Add reorder vertices.
//        2) Add ordering.hxx and add rcm and cm ordering functions.


struct Grid_connectivity_context
{
  // stupid static constexpr has external linkage.
  // while normal const does not. What is this!!!
  constexpr static int default_tag() { return 0; }

  // Observers
  // Would be required if I add operations for combinatorial editing.
  class Observer
  {
  };

  // Storage
  struct Triangle : public Removable_object
  {
    // 0 is inactive for dummies
    int vert_ids[3];
    // all active for dummies
    int tri_ids[3];
    int tag;
    // clang-format off
    Triangle(): Removable_object(), tag(default_tag())
    {vert_ids[0]=vert_ids[1]=vert_ids[2]=tri_ids[0]=tri_ids[1]=tri_ids[2]=-1;}
    // clang-format on
  };

  struct Vertex : public Removable_object
  {
    double xy[3];
    int tri_id;
    int tag;
    // clang-format off
    Vertex(): Removable_object(), tri_id(-1), tag(default_tag()) {xy[0]=xy[1]=xy[2]=0;}
    // clang-format on
  };


  Grid_connectivity_context();
  Grid_connectivity_context(const std::vector<double> & vertices,
                            const std::vector<int> & triangles,
                            const std::vector<int> & boundary_segments = std::vector<int>(),
                            const std::vector<int> & vertex_tags = std::vector<int>(),
                            const std::vector<int> & triangle_tags = std::vector<int>(),
                            const std::vector<int> & boundary_segments_tags = std::vector<int>());

  // deep copy
  Grid_connectivity_context(const Grid_connectivity_context &);

  void init(const std::vector<double> & vertices,
            const std::vector<int> & triangles,
            const std::vector<int> & boundary_segments,
            const std::vector<int> & vertex_tags,
            const std::vector<int> & triangle_tags,
            const std::vector<int> & boundary_segments_tags);
  void set_equal(const Grid_connectivity_context & g);
  void swap(Grid_connectivity_context & g);

  // Reorder vertices
  void reorder_vertices(){assert(0); }

  // Counts
  int n_vertices() const { return (int)_vertices.size(); }
  int n_real_triangles() const { return (int)(_triangles.size()) - _n_dummy_triangles; }
  int n_dummy_triangles() const { return _n_dummy_triangles; }
  int n_all_triangles() const { return (int)_triangles.size(); }

  // Access to data
  const std::vector<Vertex> & vertices_data() const { return _vertices; }
  const std::vector<Triangle> & triangles_data() const { return _triangles; }
  const std::vector<int> & boundary_tags() const { return _boundary_tags; }

private:
  static void uniqify(std::vector<int> & tags);

  int _n_dummy_triangles;
  std::vector<Vertex> _vertices;
  std::vector<Triangle> _triangles;
  std::vector<Observer *> _observers;
  std::vector<int> _boundary_tags;
};

#endif /* GRID_CONNECTIVITY_CONTEXT */
