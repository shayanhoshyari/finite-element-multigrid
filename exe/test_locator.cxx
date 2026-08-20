#include <cassert>
#include <cstdio>

#include "../src/adept/timer.h"

#include "../src/grid_connectivity_context.hxx"
#include "../src/grid_io.hxx"
#include "../src/grid_iterators.hxx"
#include "../src/point_locator.hxx"

namespace tests
{

bool
point_locator(FILE * fl)
{
  typedef bridson::Vec2d Vec2d;
  typedef bridson::Vec3d Vec3d;
  typedef bridson::Vec4d Vec4d;
  using bridson::mat_prod;

  int ierr;
  const double tol = 1e-8;

  //
  // Test triangle static function
  //
  {
    const std::vector<Vec2d> pts = {Vec2d(0, 0), Vec2d(1, 0), Vec2d(0, 1)};
    const Vec4d tr(6, 2, 7, 3);
    const Vec2d off(0.4, -1.1);

    auto run_test = [&](const Vec3d xi_exact) -> void {
      const std::vector<Vec2d> ptstr = {
        mat_prod(tr, pts[0]) + off, mat_prod(tr, pts[1]) + off, mat_prod(tr, pts[2]) + off};
      const Vec2d pos = xi_exact[0] * ptstr[0] + xi_exact[1] * ptstr[1] + xi_exact[2] * ptstr[2];
      Vec3d xi_numeric;
      ierr = Point_locator::triangle_search(pos, ptstr[0], ptstr[1], ptstr[2], xi_numeric);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      fprintf(fl,
              "e:%s, n:%s, error: %g \n",
              xi_exact.to_string().c_str(),
              xi_numeric.to_string().c_str(),
              bridson::mag(xi_exact - xi_numeric));
      assert(bridson::mag(xi_exact - xi_numeric) < tol);
    };

    run_test(Vec3d(0.3, 0.55, 0.15));
    run_test(Vec3d(-0.3, 0.55, 0.15 + 2 * 0.3));
    run_test(Vec3d(0.3, -0.55, 0.15 + 2 * 0.55));
    run_test(Vec3d(0.3, 0.55 + 2 * 0.15, -0.15));
    run_test(Vec3d(0, 1, 0));
    run_test(Vec3d(0, -1, 2));
    run_test(Vec3d(0, 0, 1));
    run_test(Vec3d(2, 0, -1));
    run_test(Vec3d(1, 0, 0));
    run_test(Vec3d(-1, 0, 2));
    run_test(Vec3d(0.5, 0.5, 0));
    fprintf(fl, "\n");
  }

  //
  // Test edge static function
  //
  {
    const std::vector<Vec2d> pts = {Vec2d(0, 0), Vec2d(1, 0)};
    // const Vec4d tr(5 * cos(.15), -sin(.15), sin(.15), 5*cos(.15));
    // const Vec2d off(0.4, -1.1);
    const Vec4d tr(4, 0, 0, 4);
    const Vec2d off(1., 1.);

    auto run_test = [&](const double x, const double y) -> void {
      const std::vector<Vec2d> ptstr = {mat_prod(tr, pts[0]) + off, mat_prod(tr, pts[1]) + off};
      const Vec2d ptrtr = mat_prod(tr, Vec2d(x, y)) + off;
      //
      const double xi_exact = std::max(0., std::min(1., x));
      const Vec2d projtr = (1 - xi_exact) * ptstr[0] + xi_exact * ptstr[1];
      const double dist_exact = bridson::mag(projtr - ptrtr);
      //
      double xi_numeric;
      double dist_numeric;
      ierr = Point_locator::edge_search(ptrtr, ptstr[0], ptstr[1], xi_numeric, dist_numeric);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      dist_numeric = sqrt(dist_numeric);
      //
      const double diff = bridson::mag(Vec2d(xi_exact - xi_numeric, dist_exact - dist_numeric));
      fprintf(fl, "e:(%g,%g), n:(%g,%g), error: %g \n", xi_exact, dist_exact, xi_numeric, dist_numeric, diff);
      // assert( diff < tol);
    };

    run_test(-1, 0);
    run_test(0, 0);

    run_test(1, 0);
    run_test(1.5, 0);
    run_test(-1, 0.33);
    run_test(0, 0.33);
    run_test(0.3, 0.33);
    run_test(1, 0.33);
    run_test(1.5, 0.33);
    fprintf(fl, "\n");
  }

  //
  // Test inside locator
  //
    {
      Point_locator::turn_on_debugging();

      Grid_connectivity_context g;
      Grid_IO(g).read_gmsh("meshes/square.msh");
      Point_locator p(g);

      Point_locator::Tri_query_in qi;
      Point_locator::Tri_query_out qo;

      qi = Point_locator::Tri_query_in(Vec2d(0.41, 0.75), 18, 6);
      ierr = p.bfs_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.directional_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.bruteforce_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);

      qi = Point_locator::Tri_query_in(Vec2d(0.41, 0.75), 86, 6);
      ierr = p.bfs_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.directional_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.bruteforce_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);

      qi = Point_locator::Tri_query_in(Vec2d(0.32, 0.92), 216, 20);
      ierr = p.bfs_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.directional_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.bruteforce_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);

      qi = Point_locator::Tri_query_in(Vec2d(0.32, 0.92), 218, 20);
      ierr = p.bfs_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.directional_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      ierr = p.bruteforce_inside_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
    }


  //
  // Test boundary locator
  //
    {
      Point_locator::turn_on_debugging();

      Grid_connectivity_context g;
      Grid_IO(g).read_gmsh("meshes/square.msh");
      Point_locator p(g);

      Point_locator::Edge_query_in qi;
      Point_locator::Edge_query_out qo;

      // No tag given
      qi = Point_locator::Edge_query_in(Vec2d(0.12, 1.11), -1, 32, 20);
      ierr = p.bfs_boundary_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      assert(qo.bedge_id == 31);
      assert(std::abs(qo.xi - .2) <= tol);
      //
      ierr = p.bruteforce_boundary_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      assert(qo.bedge_id == 31);
      assert(std::abs(qo.xi - .2) <= tol);

      // Far tag given
      qi = Point_locator::Edge_query_in(Vec2d(0.12, 1.11), 104, 32, 20);
      ierr = p.bfs_boundary_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      assert(qo.bedge_id == 7);
      assert(std::abs(qo.xi - 1) <= tol);
      //
      ierr = p.bruteforce_boundary_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      assert(qo.bedge_id == 7);
      assert(std::abs(qo.xi - 1) <= tol);

      // Near tag given
      qi = Point_locator::Edge_query_in(Vec2d(0.12, 1.11), 103, 32, 20);
      ierr = p.bfs_boundary_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      assert(qo.bedge_id == 31);
      assert(std::abs(qo.xi - .2) <= tol);
      //
      ierr = p.bruteforce_boundary_search(qi, qo);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      assert(qo.bedge_id == 31);
      assert(std::abs(qo.xi - .2) <= tol);
    }

  //
  // Test mesh locator
  //
  {

    // Locate the coarse mesh on the fine mesh
    auto test_operations = [&ierr, &tol](const Grid_connectivity_context & mesh_base,
                                         const Grid_connectivity_context & mesh_projectee,
                                         Timer & timer,
                                         const int tfast,
                                         const int texact) -> void //
    {
      Point_locator locator_base(mesh_base), locator_projectee(mesh_projectee);
      Grid_iterators giter_base(mesh_base), giter_projectee(mesh_projectee);
      Point_locator::Grid_query_out qo_exact, qo_fast;
      //
      timer.start(texact);
      ierr = locator_base.bruteforce_mesh_search(mesh_projectee, qo_exact);
      timer.stop();
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
      //
      timer.start(tfast);
      ierr = locator_base.bfs_mesh_search(mesh_projectee, qo_fast);
      timer.stop();
      assert(ierr == Point_locator::LOCATOR_SUCCESS);

      //
      // Print
      //
      const int width = 10;
      // clang-format off
      fprintf(stderr, "%*s %*s %*s %*s %*s %*s %*s\n",
                      width, "vertex",
                      width, "bdry",
                      width, "owner-bf",
                      2*width+1, "xi-bf",
                      width, "owner-fa",
                      2*width+1, "xi-fa",
                      width, "diff(fast)");
      // clang-format on

      for(int i = 0; i < giter_projectee.n_vertices(); ++i)
        {
          Grid_iterators::Vertex viter = giter_projectee.vertex_iterator(i);
          if(qo_exact.owner_id[i] == qo_fast.owner_id[i])
            {
              assert(std::abs(qo_exact.xi[2 * i + 0] - qo_fast.xi[2 * i + 0]) < tol);
              assert(std::abs(qo_exact.xi[2 * i + 1] - qo_fast.xi[2 * i + 1]) < tol);
            }
          assert(qo_exact.xi[2 * i + 0] > -tol);
          assert(qo_fast.xi[2 * i + 0] > -tol);

          double diff;
          if(!viter.is_boundary())
            {
              assert(qo_exact.xi[2 * i + 1] > -tol);
              assert(qo_fast.xi[2 * i + 1] > -tol);
              assert(1 - qo_exact.xi[2 * i + 0] - qo_exact.xi[2 * i + 1] > -tol);
              assert(1 - qo_fast.xi[2 * i + 0] - qo_fast.xi[2 * i + 1] > -tol);

              const double xi2 = 1 - qo_fast.xi[2 * i + 0] - qo_fast.xi[2 * i + 1];
              const double xi1 = qo_fast.xi[2 * i + 1];
              const double xi0 = qo_fast.xi[2 * i + 0];
              Grid_iterators::Vertex v_base0 = giter_base.triangle_iterator(qo_fast.owner_id[i]).vertex_from_offset(0);
              Grid_iterators::Vertex v_base1 = giter_base.triangle_iterator(qo_fast.owner_id[i]).vertex_from_offset(1);
              Grid_iterators::Vertex v_base2 = giter_base.triangle_iterator(qo_fast.owner_id[i]).vertex_from_offset(2);
              Grid_iterators::Vertex v_projectee = giter_projectee.vertex_iterator(i);

              const bridson::Vec2d pos_projectee(v_projectee.x(), v_projectee.y());
              const bridson::Vec2d pos_base(xi0 * v_base0.x() + xi1 * v_base1.x() + xi2 * v_base2.x(),
                                            xi0 * v_base0.y() + xi1 * v_base1.y() + xi2 * v_base2.y());
              diff = bridson::mag(pos_base - pos_projectee);
              assert(diff < tol);
            }
          else
            {
              assert(qo_exact.xi[2 * i + 0] < 1 + tol);

              const double xi1 = qo_fast.xi[2 * i + 0];
              const double xi0 = 1 - qo_fast.xi[2 * i + 0];
              Grid_iterators::Vertex v_base0 = giter_base.triangle_iterator(qo_fast.owner_id[i]).vertex_from_offset(1);
              Grid_iterators::Vertex v_base1 = giter_base.triangle_iterator(qo_fast.owner_id[i]).vertex_from_offset(2);
              Grid_iterators::Vertex v_projectee = giter_projectee.vertex_iterator(i);

              const bridson::Vec2d pos_projectee(v_projectee.x(), v_projectee.y());
              const bridson::Vec2d pos_base(xi0 * v_base0.x() + xi1 * v_base1.x(),
                                            xi0 * v_base0.y() + xi1 * v_base1.y());
              diff = bridson::mag(pos_base - pos_projectee);
              assert(diff < tol);
            }

          // clang-format off
          fprintf(stderr, "%*d %*d %*d %*g,%*g %*d %*g,%*g %*g\n",
                          width, i,
                          width, viter.is_boundary(),
                          width,  qo_exact.owner_id[i],
                          width,  qo_exact.xi[2*i+0], width,  qo_exact.xi[2*i+1],
                          width,  qo_fast.owner_id[i],
                          width,  qo_fast.xi[2*i+0], width,  qo_fast.xi[2*i+1],
                          width,   diff);
          // clang-format on

        } // End of print
    };    // End of test operations


    Point_locator::turn_off_debugging();

    Grid_connectivity_context mesh_fine, mesh_coarse, mesh_finer;
    Grid_IO(mesh_coarse).read_triangle("meshes/square.2");
    Grid_IO(mesh_fine).read_triangle("meshes/square.3");
    Grid_IO(mesh_finer).read_triangle("meshes/square.4");
    //
    Grid_IO(mesh_coarse).write_vtk("__square.2");
    Grid_IO(mesh_fine).write_vtk("__square.3");
    Grid_IO(mesh_finer).write_vtk("__square.4");

    Timer timer;
    timer.print_on_exit();

    const int t0f = timer.new_activity("COARSE TO FINE -- FAST");
    const int t0e = timer.new_activity("COARSE TO FINE -- EXACT");
    fprintf(stderr, "*** COARSE TO FINE PROJECTION STATUS \n");
    test_operations(mesh_fine, mesh_coarse, timer, t0f, t0e);

    const int t3f = timer.new_activity("FINE TO FINER -- FAST");
    const int t3e = timer.new_activity("FINE TO FINER -- EXACT");
    fprintf(stderr, "*** FINE TO FINER PROJECTION STATUS \n");
    test_operations(mesh_finer, mesh_fine, timer, t3f, t3e);

    const int t1f = timer.new_activity("FINE TO COARSE -- FAST");
    const int t1e = timer.new_activity("FINE TO COARSE -- EXACT");
    fprintf(stderr, "*** FINE TO COARSE PROJECTION STATUS \n");
    test_operations(mesh_coarse, mesh_fine, timer, t1f, t1e);

    const int t2f = timer.new_activity("FINER TO FINE -- FAST");
    const int t2e = timer.new_activity("FINER TO FINE -- EXACT");
    fprintf(stderr, "*** FINER TO FINE PROJECTION STATUS \n");
    test_operations(mesh_fine, mesh_finer, timer, t2f, t2e);


  } // End of mesh locator


  return true;
} // End of test function

} // End of namespace tests