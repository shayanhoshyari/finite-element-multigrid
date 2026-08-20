#include <algorithm>
#include <tuple>


#include "grid_connectivity_context.hxx"
#include "grid_iterators.hxx"

static int _plus_1[] = {1, 2, 0};

// ======================= Construction
Grid_connectivity_context::Grid_connectivity_context()
: _n_dummy_triangles(0)
, _vertices()
, _triangles()
, _observers()
, _boundary_tags()
{
}

Grid_connectivity_context::Grid_connectivity_context(const std::vector<double> & vertices,
                                                     const std::vector<int> & triangles,
                                                     const std::vector<int> & boundary_segments,
                                                     const std::vector<int> & vertex_tags,
                                                     const std::vector<int> & triangle_tags,
                                                     const std::vector<int> & boundary_segment_tags)
: _n_dummy_triangles(0)
, _vertices()
, _triangles()
, _observers()
, _boundary_tags()
{
  this->init(vertices, triangles, boundary_segments, vertex_tags, triangle_tags, boundary_segment_tags);
}

Grid_connectivity_context::Grid_connectivity_context(const Grid_connectivity_context & g)
: _n_dummy_triangles(0)
, _vertices()
, _triangles()
, _observers()
, _boundary_tags()
{
  this->set_equal(g);
}


void
Grid_connectivity_context::init(const std::vector<double> & vertices,
                                const std::vector<int> & triangles,
                                const std::vector<int> & boundary_segments,
                                const std::vector<int> & vertex_tags,
                                const std::vector<int> & triangle_tags,
                                const std::vector<int> & boundary_segments_tags)
{
  constexpr int tri_bs = 3; // block-size (not be mistaken with bullshit)
  constexpr int vert_bs = 3;
  constexpr int seg_bs = 2;

  assert(triangles.size() % tri_bs == 0);
  assert(vertices.size() % vert_bs == 0);
  assert(boundary_segments.size() % seg_bs == 0);

  const int num_triangles = triangles.size() / tri_bs;
  const int num_vertices = vertices.size() / vert_bs;
  const int num_segs = boundary_segments.size() / seg_bs;

  assert((vertex_tags.size() == 0) || ((int)vertex_tags.size() == num_vertices));
  assert((triangle_tags.size() == 0) || ((int)triangle_tags.size() == num_triangles));
  assert((boundary_segments_tags.size() == 0) || ((int)boundary_segments_tags.size() == num_segs));


  //
  // Add the vertices
  //
  _vertices.resize(num_vertices);
  std::fill(_vertices.begin(), _vertices.end(), Grid_connectivity_context::Vertex());
  for(int vid = 0; vid < num_vertices; vid++)
    {
      _vertices[vid].xy[0] = vertices[vid * vert_bs + 0];
      _vertices[vid].xy[1] = vertices[vid * vert_bs + 1];
      _vertices[vid].xy[2] = vertices[vid * vert_bs + 2];
      if(vertex_tags.size())
        {
          _vertices[vid].tag = vertex_tags[vid];
        }
    }

  //
  // Add the triangles
  //
  _triangles.reserve(num_triangles + num_segs);
  _triangles.resize(num_triangles);
  std::fill(_triangles.begin(), _triangles.end(), Grid_connectivity_context::Triangle());
  for(int tid = 0; tid < num_triangles; tid++)
    {
      _triangles[tid].vert_ids[0] = triangles[tid * tri_bs + 0];
      _triangles[tid].vert_ids[1] = triangles[tid * tri_bs + 1];
      _triangles[tid].vert_ids[2] = triangles[tid * tri_bs + 2];
      if(triangle_tags.size())
        {
          _triangles[tid].tag = triangle_tags[tid];
        }
      _vertices[_triangles[tid].vert_ids[0]].tri_id = tid;
      _vertices[_triangles[tid].vert_ids[1]].tri_id = tid;
      _vertices[_triangles[tid].vert_ids[2]].tri_id = tid;
    }

  //
  // Add the edges to a vector and sort them
  //
  typedef std::tuple<int, int, bool, int> Directed_edge;
  auto Directed_edge_is_boundary = [](const Directed_edge & i) -> bool { return std::get<2>(i); };
  auto Directed_edge_owner_id = [](const Directed_edge & i) -> int { return std::get<3>(i); };
  auto Directed_edge_vert_id0 = [](const Directed_edge & i) -> int { return std::get<0>(i); };
  auto Directed_edge_vert_id1 = [](const Directed_edge & i) -> int { return std::get<1>(i); };

  std::vector<Directed_edge> directed_edges;
  std::vector<int> directed_edge_sorted_ids;

  {
    directed_edges.reserve(boundary_segments.size() + num_triangles * 3);
    directed_edge_sorted_ids.reserve(boundary_segments.size() + num_triangles * 3);

    for(int tri_id = 0; tri_id < num_triangles; ++tri_id)
      {
        const int vid0 = triangles[tri_id * tri_bs + 0];
        const int vid1 = triangles[tri_id * tri_bs + 1];
        const int vid2 = triangles[tri_id * tri_bs + 2];
        directed_edges.push_back(std::make_tuple(vid0, vid1, false, tri_id));
        directed_edges.push_back(std::make_tuple(vid1, vid2, false, tri_id));
        directed_edges.push_back(std::make_tuple(vid2, vid0, false, tri_id));
        directed_edge_sorted_ids.push_back(directed_edge_sorted_ids.size());
        directed_edge_sorted_ids.push_back(directed_edge_sorted_ids.size());
        directed_edge_sorted_ids.push_back(directed_edge_sorted_ids.size());
      }
    for(int seg_id = 0; seg_id < num_segs; ++seg_id)
      {
        const int vid0 = boundary_segments[seg_id * seg_bs + 0];
        const int vid1 = boundary_segments[seg_id * seg_bs + 1];
        // I don't like this use-negative-as-segment-id-idea
        directed_edges.push_back(std::make_tuple(vid0, vid1, true, seg_id));
        directed_edge_sorted_ids.push_back(directed_edge_sorted_ids.size());
      }

    auto sort_functor = [&directed_edges](const int i, const int j) -> bool {
      const int a0 = std::get<0>(directed_edges[i]);
      const int a1 = std::get<1>(directed_edges[i]);
      const int b0 = std::get<0>(directed_edges[j]);
      const int b1 = std::get<1>(directed_edges[j]);
      const int amin = std::min(a0, a1);
      const int bmin = std::min(b0, b1);
      const int amax = std::max(a0, a1);
      const int bmax = std::max(b0, b1);

      if(amin < bmin)
        return true;
      else if(amin > bmin)
        return false;
      else if(amax < bmax)
        return true;
      else
        return false;
    };

    std::sort(directed_edge_sorted_ids.begin(), directed_edge_sorted_ids.end(), sort_functor);
  }

  //
  // Now extract the connectivity
  //
  _n_dummy_triangles = 0;

  {
    // Check if two directed edges (hyper-volumes) are equal.
    // In higher dimensions, indices should b equal in a permutation
    // invariant manner. assuming mesh is manifold, nice and ..., people
    // usually sort the indices and compare them one by one.
    auto is_equal = [](const Directed_edge & i, const Directed_edge & j) -> bool {
      const int a0 = std::get<0>(i);
      const int a1 = std::get<1>(i);
      const int b0 = std::get<0>(j);
      const int b1 = std::get<1>(j);
      const int amin = std::min(a0, a1);
      const int bmin = std::min(b0, b1);
      const int amax = std::max(a0, a1);
      const int bmax = std::max(b0, b1);
      return (amin == bmin) && (amax == bmax);
    };

    // Alternatively, I could save the offset in directed edge, i.e.,
    // typdef std::tuple<int,int,int,int> Directed_edge,
    // where entry get<3>() would be the offset.
    // For a triangular mesh this does not really matter ( I guess :) ).
    auto find_neighbour_offset = [&triangles](const int t, const int v0, const int v1) -> int {
      if((triangles[t * tri_bs + 0] == v0) && (triangles[t * tri_bs + 1] == v1)) return 2;
      if((triangles[t * tri_bs + 1] == v0) && (triangles[t * tri_bs + 0] == v1)) return 2;
      if((triangles[t * tri_bs + 1] == v0) && (triangles[t * tri_bs + 2] == v1)) return 0;
      if((triangles[t * tri_bs + 2] == v0) && (triangles[t * tri_bs + 1] == v1)) return 0;
      if((triangles[t * tri_bs + 2] == v0) && (triangles[t * tri_bs + 0] == v1)) return 1;
      if((triangles[t * tri_bs + 0] == v0) && (triangles[t * tri_bs + 2] == v1)) return 1;
      assert(0 && "Must not reach here");
      return false;
    };

    // Process a directed edge on the boundary.
    // Use a lambda expression, since this needs to be written 3 times(!)
    // in the if else clauses to come.
    auto process_directed_edge_on_boundary = [&](const Directed_edge & e0, Directed_edge const * const e1) {
      assert((!e1) || Directed_edge_is_boundary(*e1));
      const int tin_id = Directed_edge_owner_id(e0);
      const int tdum_id = (int)_triangles.size();
      const int v0 = Directed_edge_vert_id0(e0);
      const int v1 = Directed_edge_vert_id1(e0);
      const int tin_noffset = find_neighbour_offset(tin_id, v0, v1);

      _triangles.push_back(Triangle());
      _triangles.back().tri_ids[0] = tin_id;
      _triangles.back().vert_ids[2] = _triangles[tin_id].vert_ids[_plus_1[tin_noffset]];
      _triangles.back().vert_ids[1] = _triangles[tin_id].vert_ids[_plus_1[_plus_1[tin_noffset]]];
      if(e1 && boundary_segments_tags.size())
        {
          const int seg_id = Directed_edge_owner_id(*e1);
          _triangles.back().tag = boundary_segments_tags[seg_id];
        }
      ++_n_dummy_triangles;

      _triangles[tin_id].tri_ids[tin_noffset] = tdum_id;
    };

    // Process a directed edge which is internal.
    // This did not really needed to be a lambda expression,
    // but it is for consistency.
    auto process_directed_edge_internal = [&](const Directed_edge & e0, const Directed_edge & e1) {
      assert((!Directed_edge_is_boundary(e0)) && (!Directed_edge_is_boundary(e1)));
      const int t0_id = Directed_edge_owner_id(e0);
      const int t1_id = Directed_edge_owner_id(e1);
      const int v0 = Directed_edge_vert_id0(e1);
      const int v1 = Directed_edge_vert_id1(e1);
      const int t0_noffset = find_neighbour_offset(t0_id, v0, v1);
      const int t1_noffset = find_neighbour_offset(t1_id, v0, v1);
      //
      _triangles[t0_id].tri_ids[t0_noffset] = t1_id;
      _triangles[t1_id].tri_ids[t1_noffset] = t0_id;
    };


    int id = 0;
    while(id < (int)directed_edge_sorted_ids.size() - 1)
      {
        // No reference since we are swapping them
        Directed_edge e0 = directed_edges[directed_edge_sorted_ids[id]];
        Directed_edge e1 = directed_edges[directed_edge_sorted_ids[id + 1]];

        // Two consecutive equal edges
        if(is_equal(e0, e1))
          {
            if(Directed_edge_is_boundary(e0)) // put boundary at first
              std::swap(e0, e1);
            assert(!Directed_edge_is_boundary(e0));

            if(Directed_edge_is_boundary(e1)) // One is on the boundary
              process_directed_edge_on_boundary(e0, &e1);
            else // Both are internal
              process_directed_edge_internal(e0, e1);

            id += 2; // processed two
          }
        // Single edge
        else
          {
            process_directed_edge_on_boundary(e0, NULL);
            ++id; // processed one
          }
      } // End of directed edges

    // Take care of last entry
    if(id == (int)directed_edge_sorted_ids.size() - 1)
      {
        Directed_edge e0 = directed_edges[directed_edge_sorted_ids[id]];
        process_directed_edge_on_boundary(e0, NULL);
        ++id;
      }

  } // End of sort and build conn


  //
  // Last bit on Connectivity for dummy triangles, i.e., holes
  // n_triangles() now should be working
  //
  {
    int tri_id;

    //
    // First perform some checks
    //
#if((!defined(NDEBUG)) || defined(DEBUG))
    tri_id = 0;
    for(; tri_id < n_real_triangles(); ++tri_id)
      {
        const Triangle & tri = _triangles[tri_id];
        assert(tri.tri_ids[0] >= 0);
        assert(tri.tri_ids[1] >= 0);
        assert(tri.tri_ids[2] >= 0);
        assert(tri.vert_ids[0] >= 0);
        assert(tri.vert_ids[1] >= 0);
        assert(tri.vert_ids[2] >= 0);
      }
    for(; tri_id < n_all_triangles(); ++tri_id)
      {
        const Triangle & tri = _triangles[tri_id];
        assert(tri.tri_ids[0] >= 0);
        assert(tri.tri_ids[1] < 0);  // not set yet
        assert(tri.tri_ids[2] < 0);  // not set yet
        assert(tri.vert_ids[0] < 0); // Will never be set
        assert(tri.vert_ids[1] >= 0);
        assert(tri.vert_ids[2] >= 0);
      }
#endif

    //
    // Now set the last bit
    //
    Grid_iterators grid_iterators(*this);

    tri_id = n_real_triangles();
    for(; tri_id < n_all_triangles(); tri_id++)
      {
        Triangle & tri_data = _triangles[tri_id];
        Grid_iterators::Half_edge he = grid_iterators.half_edge_iterator(tri_id, 0);
        Grid_iterators::Half_edge he_next = he;

        //   Set the next
        he_next = he.twin();
        while(!he_next.is_boundary())
          {
            he_next = he_next.next().next().twin();
          }
        tri_data.tri_ids[1] = he_next.triangle().raw_id();

        // Set the prev
        he_next = he.twin();
        while(!he_next.is_boundary())
          {
            he_next = he_next.next().twin();
          }
        tri_data.tri_ids[2] = he_next.triangle().raw_id();

      } // End of hald edges
  }     // End of fixing holes


  //
  // A last check that all is well
  //
  {
#if((!defined(NDEBUG)) || defined(DEBUG))
    int tri_id = n_real_triangles();
    for(; tri_id < n_all_triangles(); ++tri_id)
      {
        const Triangle & tri = _triangles[tri_id];
        assert((tri.tri_ids[0] >= 0) && (tri.tri_ids[0] < n_real_triangles()));
        assert((tri.tri_ids[1] >= n_real_triangles()) && (tri.tri_ids[1] != tri_id)); // must now be set
        assert((tri.tri_ids[2] >= n_real_triangles()) && (tri.tri_ids[2] != tri_id)); // must now be set
        assert(tri.vert_ids[0] < 0);                                                  // Will never be set
        assert(tri.vert_ids[1] >= 0);
        assert(tri.vert_ids[2] >= 0);
      }
#endif
  } // End of final check


  // Now get the boundary tags
  if(boundary_segments_tags.size())
    {
      _boundary_tags = boundary_segments_tags;
      Grid_connectivity_context::uniqify(_boundary_tags);
    }
  else
    {
      _boundary_tags = std::vector<int>(1, default_tag());
    }

} // All done

void
Grid_connectivity_context::set_equal(const Grid_connectivity_context & g)
{
  _n_dummy_triangles = g._n_dummy_triangles;
  _vertices = g._vertices;
  _triangles = g._triangles;
  _boundary_tags = g._boundary_tags;
  // I don't care about observers now.
}

void
Grid_connectivity_context::swap(Grid_connectivity_context & g)
{
  _n_dummy_triangles = g._n_dummy_triangles;
  _vertices.swap(g._vertices);
  _triangles.swap(g._triangles);
  _boundary_tags.swap(g._boundary_tags);
  // I don't care about observers now.
}

void
Grid_connectivity_context::uniqify(std::vector<int> & tags)
{
  std::sort(tags.begin(), tags.end());
  std::vector<int>::iterator it_last = std::unique(tags.begin(), tags.end());
  tags.erase(it_last, tags.end());
}
