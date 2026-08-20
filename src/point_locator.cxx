#include <cmath>
#include <queue>
#include <set>
#include <tuple>

#include "bridson/vec.h"

#include "fake_fprintf.hxx"
#include "grid_iterators.hxx"
#include "point_locator.hxx"

// clang-format off
// #define CHKERRQ(ierr) assert(ierr != LOCATOR_FATAL_ERROR)
#define CHKERRQ(ierr)  do{ if(ierr == LOCATOR_FATAL_ERROR) {return LOCATOR_FATAL_ERROR;}}while(0)
#define FL stderr
#if !defined(NDEBUG)
#   define FPRINTF   Point_locator::_printf_functor
#else
#   define FPRINTF   fake_fprintf
#endif
// clang-format on

Fprintf_function Point_locator::_printf_functor = fake_fprintf;
void
Point_locator::turn_on_debugging()
{
  _printf_functor = fprintf;
}
void
Point_locator::turn_off_debugging()
{
  _printf_functor = fake_fprintf;
}


Point_locator::Point_locator(const Grid_connectivity_context & conn)
: _connectivity(&conn)
{
}

static double
cross2d(const bridson::Vec2d & a, const bridson::Vec2d & b)
{
  return a[0] * b[1] - a[1] * b[0];
}

// ==============================================================
//            Search a signle triangle
// ==============================================================
int
Point_locator::triangle_search(const Vec2d & x, const Vec2d & p0, const Vec2d & p1, const Vec2d & p2, Vec3d & xi)
{
  const bridson::Vec2d e1 = p1 - p0;
  const bridson::Vec2d e2 = p2 - p0;
  const bridson::Vec2d e3 = p0 - x;
  const bridson::Vec2d e4 = p1 - x;
  const bridson::Vec2d e5 = p2 - x;

  const double Atot = cross2d(e1, e2);
  const double A0 = cross2d(e4, e5);
  const double A1 = cross2d(e5, e3);
  const double A2 = cross2d(e3, e4);

  xi[0] = A0 / Atot;
  xi[1] = A1 / Atot;
  xi[2] = A2 / Atot;

  if(std::abs(1. - xi[0] - xi[1] - xi[2]) > tol())
    {
      return LOCATOR_FATAL_ERROR;
    }
  return LOCATOR_SUCCESS;
}

int
Point_locator::triangle_search(const Vec2d & x, const int tri_id, Vec3d & xi) const
{
  const Grid_iterators giter(connectivity());
  const Grid_iterators::Triangle tri = giter.triangle_iterator(tri_id);
  assert(!tri.is_dummy());
  const Grid_iterators::Vertex v0 = tri.vertex_from_offset(0);
  const Grid_iterators::Vertex v1 = tri.vertex_from_offset(1);
  const Grid_iterators::Vertex v2 = tri.vertex_from_offset(2);
  const Vec2d p0(v0.x(), v0.y());
  const Vec2d p1(v1.x(), v1.y());
  const Vec2d p2(v2.x(), v2.y());
  return Point_locator::triangle_search(x, p0, p1, p2, xi);
}

// ==============================================================
//            Project On a single edge
// ==============================================================

int
Point_locator::edge_search(const Vec2d & x, const Vec2d & p0, const Vec2d & p1, double & xi, double & dist2)
{
  const double tol = 1e-11;

  const Vec2d e1 = (p1 - p0);
  const Vec2d e2 = (x - p0);
  const double e1_mag2 = bridson::mag2(e1);
  {
    assert(e1_mag2 > tol);
    const double tt = bridson::dot(e1, e2) / e1_mag2;
    xi = std::max(std::min(tt, 1.), 0.);
  }
  const Vec2d closest_point = (1 - xi) * p0 + xi * p1;
  dist2 = bridson::mag2(closest_point - x);
  return LOCATOR_SUCCESS;
}

int
Point_locator::boundary_edge_search(const Vec2d & x, const int bedge_id, double & xi, double & dist2) const
{
  const Grid_iterators giter(connectivity());
  const Grid_iterators::Half_edge he = giter.boundary_half_edge_iterator(bedge_id);
  const Grid_iterators::Vertex v0 = he.origin();
  const Grid_iterators::Vertex v1 = he.twin().origin();
  const Vec2d p0(v0.x(), v0.y());
  const Vec2d p1(v1.x(), v1.y());
  return Point_locator::edge_search(x, p0, p1, xi, dist2);
}


// ==============================================================
//           Querry constructors
// ==============================================================

Point_locator::Tri_query_in::Tri_query_in(const Vec2d _x, const int _tri_id, const int _max_bfs)
: x(_x)
, tri_id(_tri_id)
, max_bfs_layers(_max_bfs)
{
}

Point_locator::Tri_query_out::Tri_query_out(const Vec3d _xi, const int _tri_id)
: xi(_xi)
, tri_id(_tri_id)
{
}

Point_locator::Edge_query_in::Edge_query_in(const Vec2d _x, const int _tag, const int _edge, const int _max_bfs)
: x(_x)
, tag(_tag)
, bedge_id(_edge)
, max_bfs_layers(_max_bfs)
{
}

Point_locator::Edge_query_out::Edge_query_out(const double _xi, const int _edge, const double _dist2)
: xi(_xi)
, bedge_id(_edge)
, dist2(_dist2)
{
}

// ==============================================================
//            Find a point inside a mesh
// ==============================================================

//
// Partially tested. Much faster than BFS.
//
int
Point_locator::directional_inside_search(Tri_query_in i, Tri_query_out & o) const
{
  int n_attempts = 0;
  std::set<int> visited;
  FPRINTF(FL, "* BEGIN DFS INTERNAL search begin, pt:%s guess:%d\n", i.x.to_string().c_str(), i.tri_id);
  int ierr = this->directional_inside_search_helper(i, o, n_attempts, visited);
  FPRINTF(FL, "....dfs directional search ended ierr=%d, n_attempts=%d\n", ierr, n_attempts);
  return ierr;
}

//
// Usually the record of visited people is not needed.
// But it is added to prevent pathological unending searches.
//
int
Point_locator::directional_inside_search_helper(Tri_query_in i,
                                                Tri_query_out & o,
                                                int & n_attempts,
                                                std::set<int> & visited) const
{
  Grid_iterators g(connectivity());
  int ierr;
  ++n_attempts;

  Grid_iterators::Triangle tri_guess = g.triangle_iterator(i.tri_id);

  assert(!visited.count(i.tri_id));
  visited.insert(i.tri_id);
  ierr = triangle_search(i.x, i.tri_id, o.xi);
  CHKERRQ(ierr);
  if(bridson::min(o.xi) > -tol())
    {
      o.tri_id = i.tri_id;
      FPRINTF(FL, "....dfs internal search found at tri=%d, xi=%s \n", tri_guess.raw_id(), o.xi.to_string().c_str());
      return LOCATOR_SUCCESS;
    }


  int next_faces[] = {tri_guess.neighbour_from_offset(0).raw_id(),
                      tri_guess.neighbour_from_offset(1).raw_id(),
                      tri_guess.neighbour_from_offset(2).raw_id()};
  int is_in[] = {!tri_guess.neighbour_from_offset(0).is_dummy(),
                 !tri_guess.neighbour_from_offset(1).is_dummy(),
                 !tri_guess.neighbour_from_offset(2).is_dummy()};

  if((o.xi[0] < 0) && (is_in[0]) && (!visited.count(next_faces[0])))
    {
      i.tri_id = next_faces[0];
      FPRINTF(FL, "....dfs internal  search xi0<0 visiting %d \n", i.tri_id);
      return this->directional_inside_search_helper(i, o, n_attempts, visited);
    }
  else if((o.xi[1] < 0) && (is_in[1]) && (!visited.count(next_faces[1])))
    {
      i.tri_id = next_faces[1];
      FPRINTF(FL, "....dfs internal  search xi1<0 visiting %d \n", i.tri_id);
      return this->directional_inside_search_helper(i, o, n_attempts, visited);
    }
  else if((o.xi[2] < 0) && (is_in[2]) && (!visited.count(next_faces[2])))
    {
      i.tri_id = next_faces[2];
      FPRINTF(FL, "....dfs internal  search xi2<0 visiting %d \n", i.tri_id);
      return this->directional_inside_search_helper(i, o, n_attempts, visited);
    }

  return LOCATOR_FAILURE;
}

//
// Partially tested.
//
int
Point_locator::bfs_inside_search(Tri_query_in ii, Tri_query_out & oo) const
{
  int ierr;
  int n_attempts = 0;
  Grid_iterators giter(connectivity());

  typedef std::pair<int, int> Index_and_depth;
  auto get_index = [](const Index_and_depth & a) -> int { return a.first; };
  auto get_depth = [](const Index_and_depth & a) -> int { return a.second; };

  std::queue<Index_and_depth> to_inspect;
  // If the number of layers is small, can also use a vector
  std::set<int> have_visited;

  const int max_depth = ii.max_bfs_layers;
  to_inspect.push(std::make_pair(ii.tri_id, 0));
  assert(!giter.triangle_iterator(ii.tri_id).is_dummy());

  FPRINTF(FL, "* BEGIN BFS INTERNAL SEARCH Point: %s guess: %d \n", ii.x.to_string().c_str(), ii.tri_id);

  //
  // Start the BFS search
  //
  while(!to_inspect.empty())
    {
      // Pop the element
      Index_and_depth index_and_depth = to_inspect.front();
      to_inspect.pop();
      const int depth = get_depth(index_and_depth);
      Grid_iterators::Triangle tri = giter.triangle_iterator(get_index(index_and_depth));
      FPRINTF(FL, "....bfs_inside_search: popped %d with depth %d \n", tri.raw_id(), depth);

      // remember that it was visited
      if(have_visited.count(tri.raw_id()))
        {
          FPRINTF(FL, "... have visitied, skip \n");
          continue;
        }
      else
        {
          have_visited.insert(tri.raw_id());
        }

      // see if triangle is inside. if inside, happily return
      ierr = triangle_search(ii.x, tri.raw_id(), oo.xi);
      CHKERRQ(ierr);
      ++n_attempts;
      if(bridson::min(oo.xi) > -tol())
        {
          oo.tri_id = tri.raw_id();
          FPRINTF(FL,
                  "....bfs_inside_search: Found, tri:%d pos:%s attempts:%d \n",
                  tri.raw_id(),
                  oo.xi.to_string().c_str(),
                  n_attempts);
          return LOCATOR_SUCCESS;
        }
      FPRINTF(FL, "....bfs_inside_search: Not here, tri:%d pos:%s \n", tri.raw_id(), oo.xi.to_string().c_str());

      // it was not inside, push in neigbours
      if(depth < max_depth)
        {
          for(int offset = 0; offset < 3; ++offset)
            {
              Grid_iterators::Triangle nei = tri.neighbour_from_offset(offset);
              FPRINTF(FL, "....Neighbour %d is %d, ", offset, nei.raw_id());
              if(nei.is_dummy())
                {
                  FPRINTF(FL, "..Is boundary skip... \n");
                }
              else
                {
                  FPRINTF(FL, "....Push to queueueueue \n");
                  to_inspect.push(std::make_pair(nei.raw_id(), depth + 1));
                }
            }
        }
      else
        {
          FPRINTF(FL, "....Reached max depth. Won't check neighbours \n");
        }
    } // End of BFS search

  return LOCATOR_FAILURE;
} // All done

//
// Partially tested.
//
int
Point_locator::bruteforce_inside_search(Tri_query_in ii, Tri_query_out & oo) const
{
  int ierr;
  Grid_iterators giter(connectivity());

  FPRINTF(FL, "* BEGIN BRUTE FORCE INTERNAL SEARCH Point: %s guess: %d \n", ii.x.to_string().c_str(), ii.tri_id);

  //
  // Brute force search
  //
  for(int tri_id = 0; tri_id < giter.n_real_triangles(); ++tri_id)
    {
      Grid_iterators::Triangle tri = giter.triangle_iterator(tri_id);

      // see if triangle is inside. if inside, happily return
      ierr = triangle_search(ii.x, tri.raw_id(), oo.xi);
      CHKERRQ(ierr);
      if(bridson::min(oo.xi) > -tol())
        {
          FPRINTF(FL, "....brute force internal succeeded xi: %s tri: %d \n", oo.xi.to_string().c_str(), oo.tri_id);
          oo.tri_id = tri.raw_id();
          return LOCATOR_SUCCESS;
        }
    } // End of Brute force search

  // Brute force cannot fail, so if it does not find something, it is a fatal error.
  return LOCATOR_FATAL_ERROR;
} // All done

// ==============================================================
//            Project on the boundary of the mesh
// ==============================================================
int
Point_locator::bfs_boundary_search(Edge_query_in ii, Edge_query_out & oo) const
{
  int ierr;
  Grid_iterators giter(connectivity());
  //
  typedef std::pair<int, int> Index_and_depth;
  auto get_index = [](const Index_and_depth & a) -> int { return a.first; };
  auto get_depth = [](const Index_and_depth & a) -> int { return a.second; };
  //
  const int max_depth = ii.max_bfs_layers;
  std::queue<Index_and_depth> to_inspect;
  std::set<int> have_visited;

  FPRINTF(
    FL, "* BEGIN BFS BOUNDARY SEARCH Point: %s guess: %d tag: %d \n", ii.x.to_string().c_str(), ii.bedge_id, ii.tag);

  //  Make sure nonesense is not given to us
  assert((giter.boundary_half_edge_iterator(ii.bedge_id).specific_half_edge_id() >= 0) &&
         "boundary edge id is invalid");
  // make sure the tag is write
  assert((ii.tag < 0) || is_valid_btag(ii.tag));

  // Push the first guy into the queue
  to_inspect.push(std::make_pair(ii.bedge_id, 0));

  // Record the closest instance.
  double xi;
  double dist2;
  oo.dist2 = std::numeric_limits<double>::max();
  oo.bedge_id = -1;

  //
  // Start the BFS search
  //
  while(!to_inspect.empty())
    {
      // get the element
      Index_and_depth index_and_depth = to_inspect.front();
      to_inspect.pop();
      const int depth = get_depth(index_and_depth);
      Grid_iterators::Half_edge he = giter.boundary_half_edge_iterator(get_index(index_and_depth));
      FPRINTF(FL,
              "....bfs_boundary_search: popped bhe %d tag: %d  depth: %d ;",
              he.specific_half_edge_id(),
              he.triangle().tag(),
              depth);

      // remember and check if that it was visited
      if(have_visited.count(he.specific_half_edge_id()))
        {
          FPRINTF(FL, " have visitied, skip ...\n");
          continue;
        }
      else
        {
          have_visited.insert(he.specific_half_edge_id());
        }

      // Check the tag and then do the projection
      // But don't skip, since we want to go through the neighbours
      if((ii.tag <= 0) || (he.triangle().tag() == ii.tag))
        {
          ierr = boundary_edge_search(ii.x, he.specific_half_edge_id(), xi, dist2);
          CHKERRQ(ierr);
          FPRINTF(FL, " dist2 is %g; ", dist2);

          if(dist2 < oo.dist2)
            {
              oo.dist2 = dist2;
              oo.xi = xi;
              oo.bedge_id = he.specific_half_edge_id();
              FPRINTF(FL, "best so far");
            }
          // If point lies on boundary exactly, just get out
          if(dist2 < tol())
            {
              FPRINTF(FL, "; too good, let's get out of here \n");
              break;
            }
          FPRINTF(FL, "\n");
        }
      else
        {
          FPRINTF(FL, " tag incompatible, ignore \n");
        }

      // it was not inside, push in neigbours
      // If depth is reached, but we have found nothing, still continue to look.
      if((depth < max_depth) || (oo.bedge_id < 0))
        {
          Grid_iterators::Half_edge nei = he.next();
          FPRINTF(FL, "... neighour 0 is %d, push to queue\n", nei.specific_half_edge_id());
          to_inspect.push(std::make_pair(nei.specific_half_edge_id(), depth + 1));
          //
          nei = he.prev();
          FPRINTF(FL, "... neighour 1 is %d, push to queue\n", nei.specific_half_edge_id());
          to_inspect.push(std::make_pair(nei.specific_half_edge_id(), depth + 1));
        }
    } // End of BFS search

  // IF a wrong tag is given then fatal error occurs
  if(oo.bedge_id >= 0)
    return LOCATOR_SUCCESS;
  else
    return LOCATOR_FATAL_ERROR;
}

int
Point_locator::bruteforce_boundary_search(Edge_query_in ii, Edge_query_out & oo) const
{
  int ierr;
  Grid_iterators giter(connectivity());

  FPRINTF(FL, "* BEGIN BRUTEFORCE BOUNDARY SEARCH Point: %s guess: %d \n", ii.x.to_string().c_str(), ii.bedge_id);

  double xi;
  double dist2;
  oo.dist2 = std::numeric_limits<double>::max();
  oo.bedge_id = -1;

  //  Make sure nonesense is not given to us
  assert((giter.boundary_half_edge_iterator(ii.bedge_id).specific_half_edge_id() > 0) && "boundary edge id is invalid");
  // make sure the tag is write
  assert((ii.tag < 0) || is_valid_btag(ii.tag));

  //
  // Brute force search
  //
  for(int bedge_id = 0; bedge_id < giter.n_boundary_edges(); ++bedge_id)
    {
      // get the element
      Grid_iterators::Half_edge he = giter.boundary_half_edge_iterator(bedge_id);

      // Don't consider edges with a different tag
      // If a positive number for the tag is specified
      if((ii.tag >= 0) && he.triangle().tag() != ii.tag) continue;

      // see if triangle is inside. if inside, happily return
      ierr = boundary_edge_search(ii.x, he.specific_half_edge_id(), xi, dist2);
      CHKERRQ(ierr);
      if(dist2 < oo.dist2)
        {
          oo.dist2 = dist2;
          oo.bedge_id = he.specific_half_edge_id();
          oo.xi = xi;
        }
      // IF too good, just break;
      if(dist2 < tol())
        {
          break;
        }

    } // End of Brute force search

  // IF a wrong tag is given then fatal error occurs
  if(oo.bedge_id >= 0)
    {
      FPRINTF(FL, "....brute force boundary succeeded xi:%g edge:%d dist2:%g\n", oo.xi, oo.bedge_id, oo.dist2);
      return LOCATOR_SUCCESS;
    }
  else
    {
      FPRINTF(FL, "....brute force boundary failed \n");
      return LOCATOR_FATAL_ERROR;
    }
} // All done

// ==============================================================
//            Project one mesh on another
// ==============================================================

// Static helpers (seemed cleaner in this occasion compared to lambda's
namespace
{

struct Searchee
{
  int invertex;
  int mytri; // must be ghost for boundary
  int btag;
  Searchee(const int invertex_, const int mytri_, const int btag_)
  : invertex(invertex_)
  , mytri(mytri_)
  , btag(btag_)
  {
  }
};

static void
bfs_mesh_push_boundary_neighbours_to_queue(const Grid_iterators & in,
                                           const int v,
                                           const Point_locator::Grid_query_out & oo,
                                           std::queue<Searchee> & queue)
{

  Grid_iterators::Vertex_umbrella umbrella = in.vertex_umbrella_iterator(v);
  const int v_owner = oo.owner_id[v];
  assert(v_owner >= 0);
  //
  bool is_on_boundary;
  is_on_boundary = umbrella.reset_boundary();
  assert(is_on_boundary);

  //
  // Good old pedanctic enumeration
  //

  // Half-edges
  const Grid_iterators::Half_edge hea = umbrella.half_edge();
  const Grid_iterators::Half_edge heb = hea.next();
  const Grid_iterators::Half_edge hec = heb.twin();
  const Grid_iterators::Half_edge hed = hea.twin();
  const Grid_iterators::Half_edge hee = hea.prev();
  const Grid_iterators::Half_edge hef = heb.next();

  // Triangles
  const Grid_iterators::Triangle tria = hea.triangle();
  const Grid_iterators::Triangle trib = heb.triangle();
  // const Grid_iterators::Triangle tric = hec.triangle();
  // const Grid_iterators::Triangle trid = hed.triangle();
  const Grid_iterators::Triangle trie = hee.triangle();
  const Grid_iterators::Triangle trif = hef.triangle();

  // Vertices
  const Grid_iterators::Vertex va = hea.origin();
  const Grid_iterators::Vertex vb = hec.origin();

  assert(hea.is_boundary());
  assert(heb.is_boundary());
  assert(hee.is_boundary());
  assert(hef.is_boundary());
  //
  assert(tria.is_dummy());
  assert(trib.is_dummy());
  assert(trie.is_dummy());
  assert(trif.is_dummy());

  {
    // consistency
    const int btag = std::min(trie.tag(), tria.tag());
    queue.push(Searchee(va.id(), v_owner, btag));
  }
  {
    // consistency
    const int btag = std::min(trib.tag(), trif.tag());
    queue.push(Searchee(vb.id(), v_owner, btag));
  }
}

static void
bfs_mesh_push_internal_neighbours_to_queue(const Grid_iterators & in,
                                           const int v,
                                           const Point_locator::Grid_query_out & oo,
                                           std::queue<Searchee> & queue)
{
  Grid_iterators::Vertex_umbrella umbrella = in.vertex_umbrella_iterator(v);
  const int v_owner = oo.owner_id[v];
  assert(v_owner >= 0);
  assert(!in.vertex_iterator(v).is_boundary());
  do
    {
      Grid_iterators::Vertex nei = umbrella.half_edge().origin();
      if(!nei.is_boundary())
        {
          queue.push(Searchee(nei.id(), v_owner, -1));
        }
    }
  while(umbrella.advance());
}


} // End of anonymus napspace


void
Point_locator::bfs_mesh_boundary_propagative_search(const Grid_connectivity_context & inconn,
                                                    const int seed_vertex_id,
                                                    Point_locator::Grid_query_out & oo) const
{
  const Grid_iterators me = Grid_iterators(connectivity());
  const Grid_iterators in = Grid_iterators(inconn);
  int ierr;

  std::queue<Searchee> propagation_front;

  const int seed_owner = oo.owner_id[seed_vertex_id];
  assert(seed_owner >= 0);
  bfs_mesh_push_boundary_neighbours_to_queue(in, seed_vertex_id, oo, propagation_front);

  while(!propagation_front.empty())
    {
      Searchee searchee = propagation_front.front();
      propagation_front.pop();
      if(oo.owner_id[searchee.invertex] >= 0) continue;

      Grid_iterators::Vertex searche_v = in.vertex_iterator(searchee.invertex);
      Grid_iterators::Triangle searche_trime = me.triangle_iterator(searchee.mytri);
      assert(searche_trime.is_dummy());

      Edge_query_in qi(bridson::Vec2d(searche_v.x(), searche_v.y()),
                       searchee.btag,
                       searche_trime.specific_id(),
                       MAX_BOUNDARY_BFS_LAYERS);
      Edge_query_out qo;
      ierr = this->bfs_boundary_search(qi, qo);
      assert(ierr == LOCATOR_SUCCESS);

      oo.owner_id[searche_v.id()] = me.boundary_half_edge_iterator(qo.bedge_id).triangle().raw_id();
      oo.xi[2 * searche_v.id() + 0] = qo.xi;
      oo.xi[2 * searche_v.id() + 1] = -1000;

      bfs_mesh_push_boundary_neighbours_to_queue(in, searche_v.id(), oo, propagation_front);
    }
}

void
Point_locator::bfs_mesh_internal_propagative_search(const Grid_connectivity_context & inconn,
                                                    const int seed_vertex_id,
                                                    const int seed_mytri_guess,
                                                    Point_locator::Grid_query_out & oo) const
{
  const Grid_iterators me = Grid_iterators(connectivity());
  const Grid_iterators in = Grid_iterators(inconn);
  int ierr;

  std::queue<Searchee> propagation_front;
  propagation_front.push(Searchee(seed_vertex_id, seed_mytri_guess, -1));

  while(!propagation_front.empty())
    {
      Searchee searchee = propagation_front.front();
      propagation_front.pop();
      if(oo.owner_id[searchee.invertex] >= 0) continue;

      Grid_iterators::Vertex searche_v = in.vertex_iterator(searchee.invertex);
      Grid_iterators::Triangle searche_trime = me.triangle_iterator(searchee.mytri);
      assert(!searche_trime.is_dummy());

      Tri_query_in qi(bridson::Vec2d(searche_v.x(), searche_v.y()), searche_trime.raw_id(), -100 /* does not matter */);
      Tri_query_out qo;
      // Directional is faster (should I make it DFS to be more robust?)
      // qi.max_bfs_layers = MAX_BFS_LAYERS;
      // ierr = this->bfs_inside_search(qi, qo);
      ierr = this->directional_inside_search(qi, qo);
      assert(ierr == LOCATOR_SUCCESS);
      assert(!me.triangle_iterator(qo.tri_id).is_dummy());

      oo.owner_id[searche_v.id()] = qo.tri_id;
      oo.xi[2 * searche_v.id() + 0] = qo.xi[0];
      oo.xi[2 * searche_v.id() + 1] = qo.xi[1];

      bfs_mesh_push_internal_neighbours_to_queue(in, searche_v.id(), oo, propagation_front);
    }
}

int
Point_locator::bfs_mesh_search(const Grid_connectivity_context & inconn, Grid_query_out & oo) const
{

  Grid_iterators my(connectivity());
  Grid_iterators in(inconn);
  int ierr;

  oo.owner_id.resize(in.n_vertices());
  oo.xi.resize(in.n_vertices() * 2);
  std::fill(oo.xi.begin(), oo.xi.end(), -10000);
  std::fill(oo.owner_id.begin(), oo.owner_id.end(), -1);

  for(int vinid = 0; vinid < in.n_vertices(); ++vinid)
    {
      // Querry if we are on the boundary
      Grid_iterators::Vertex viniter = in.vertex_iterator(vinid);
      int btag0, btag1;
      const bool isbdry = viniter.is_boundary(btag0, btag1);

      // If we were indeed on the boundary and our owner is unknown, then
      // find it by brute force, and then propagate the information
      if((oo.owner_id[vinid] < 0) && isbdry)
        {
          // consistency
          const int btag = std::min(btag0, btag1);

          // Find by brute force
          Edge_query_in qi(Vec2d(viniter.x(), viniter.y()), btag, -1, -1);
          Edge_query_out qo;
          ierr = bruteforce_boundary_search(qi, qo);
          assert(ierr == LOCATOR_SUCCESS);
          Grid_iterators::Half_edge my_boundary_he_iter = my.boundary_half_edge_iterator(qo.bedge_id);
          assert(my_boundary_he_iter.is_boundary());
          oo.owner_id[viniter.id()] = my_boundary_he_iter.triangle().raw_id();
          oo.xi[2 * viniter.id() + 0] = qo.xi;
          oo.xi[2 * viniter.id() + 1] = -1000;

          // Propagate
          bfs_mesh_boundary_propagative_search(inconn, viniter.id(), oo);
        } // End of no owner && boundary

      //
      // If we are on the boundary
      // 1) assert that now we have an owner
      // 2) see if we have an internal neighour with no owner and propagate
      //
      if(isbdry)
        {
          assert(oo.owner_id[vinid] >= 0);
          const Grid_iterators::Half_edge my_boundary_he_iter = my.triangle_iterator(oo.owner_id[vinid]).half_edge();

          // Find a neighbour vertex which does not touch the boundary
          bool has_inner_neighbour = false;
          Grid_iterators::Vertex_umbrella vinumber = in.vertex_umbrella_iterator(viniter.id());
          vinumber.reset_boundary();
          do
            {
              Grid_iterators::Vertex vother = vinumber.half_edge().origin();
              if((oo.owner_id[vother.id()] < 0) && (!vother.is_boundary()))
                {
                  has_inner_neighbour = true;
                  break;
                }
            }
          while(vinumber.advance());

          // If it exists, then propagate
          if(has_inner_neighbour)
            {
              Grid_iterators::Vertex viter_inmesh_inside = vinumber.half_edge().origin();
              bfs_mesh_internal_propagative_search(
                inconn, viter_inmesh_inside.id(), my_boundary_he_iter.twin().triangle().raw_id(), oo);
            }
        } // End of is_boundary


    } // End of loop


  //
  // MAke sure all are visited
  //
  for(int vinid = 0; vinid < in.n_vertices(); ++vinid)
    {
      assert(oo.owner_id[vinid] >= 0);
    }


  return LOCATOR_SUCCESS;
}

int
Point_locator::bruteforce_mesh_search(const Grid_connectivity_context & in, Grid_query_out & oo) const
{
  int ierr;
  Grid_iterators myiters(connectivity());
  Grid_iterators initers(in);

  assert(in.boundary_tags() == connectivity().boundary_tags());

  oo.owner_id.resize(initers.n_vertices());
  oo.xi.resize(2 * initers.n_vertices());

  for(int i = 0; i < initers.n_vertices(); ++i)
    {
      Grid_iterators::Vertex vin = initers.vertex_iterator(i);
      int btag0, btag1;

      if(vin.is_boundary(btag0, btag1))
        {
          // use the minimum of btags for consistency
          const int btag = std::min(btag0, btag1);

          Edge_query_in qin(Vec2d(vin.x(), vin.y()), btag, -1, -1);
          Edge_query_out qo;
          ierr = bruteforce_boundary_search(qin, qo);
          assert(ierr == LOCATOR_SUCCESS);
          CHKERRQ(ierr);
          Grid_iterators::Half_edge he = myiters.boundary_half_edge_iterator(qo.bedge_id);
          oo.owner_id[vin.id()] = he.triangle().raw_id();
          oo.xi[2 * vin.id() + 0] = qo.xi;
          oo.xi[2 * vin.id() + 1] = -1000;
        }
      else
        {
          Tri_query_in qin(Vec2d(vin.x(), vin.y()), -1, -1);
          Tri_query_out qo;
          ierr = bruteforce_inside_search(qin, qo);
          assert(ierr == LOCATOR_SUCCESS);
          CHKERRQ(ierr);
          oo.owner_id[vin.id()] = qo.tri_id;
          oo.xi[2 * vin.id() + 0] = qo.xi[0];
          oo.xi[2 * vin.id() + 1] = qo.xi[1];
        }
    }

  return LOCATOR_SUCCESS;
}
