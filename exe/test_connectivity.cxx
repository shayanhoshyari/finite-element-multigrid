#include <cstdio>
#include <cassert>

#include "../src/grid_io.hxx"
#include "../src/grid_connectivity_context.hxx"

namespace tests
{

bool connectivity(FILE *fl)
{
    // clang-format off
  std::vector<double> pos =
  {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0
  };
  std::vector<int> conn =
  {
    0, 1, 3,
    1, 2, 3
  };
  std::vector<int> segs =
  {
    0,1,
    1,2,
    3,2,
  };
  std::vector<int> seg_tags =
  {
    21,
    23,
    33
  };
  // clang-format on

  fprintf(fl, "== TESTING CONNECTIVITY\n");

  Grid_connectivity_context g(pos, conn);
  Grid_IO(g).print_data_info(fl);
  Grid_IO(g).print_iter_info(fl);
  Grid_IO(g).write_vtk("__test.vtk");

  g.init(pos, conn, segs, std::vector<int>(), std::vector<int>(), seg_tags);
  Grid_IO(g).print_data_info(fl);
  Grid_IO(g).print_iter_info(fl);
  Grid_IO(g).write_vtk("__test2.vtk");

  //Grid_IO(g).read_obj("meshes/cartel/cow1.obj");
  //Grid_IO(g).write_vtk("__cow.vtk");

  Grid_IO(g).read_triangle("meshes/square.2");
  Grid_IO(g).write_vtk("__square.2.vtk");

  return true;
}

};