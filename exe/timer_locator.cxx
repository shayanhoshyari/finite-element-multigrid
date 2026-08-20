#include <cassert>
#include <cstdio>

#include "../src/adept/timer.h"

#include "../src/grid_connectivity_context.hxx"
#include "../src/grid_io.hxx"
#include "../src/grid_iterators.hxx"
#include "../src/point_locator.hxx"


//
// Test mesh locator
//
int
main(int argc, char * argv[])
{

  int ierr;

  int perform_exact = 1;
  std::string s0("meshes/square.4"), s1("meshes/square.5");
  if(argc >= 2) s0 = argv[1];
  if(argc >= 3) s1 = argv[2];
  if(argc >= 4) sscanf(argv[3], "%d", &perform_exact);

  printf("Mesh 1: %s\n", s0.c_str());
  printf("Mesh 2: %s\n", s1.c_str());


  Grid_connectivity_context c0, c1;
  Grid_IO io0(c0), io1(c1);

  io0.read_triangle(s0);
  io1.read_triangle(s1);
  io0.write_vtk("__proj1.vtk");
  io1.write_vtk("__proj2.vtk");

  Point_locator l0(c0), l1(c1);
  Point_locator::Grid_query_out querry;

  Timer timer;
  timer.print_on_exit();


  const int t0f = timer.new_activity("FIRST TO SECOND -- FAST");
  const int t0e = timer.new_activity("FIRST TO SECOND -- EXACT");
  fprintf(stderr, "*** FIRST TO SECOND PROJECTION STATUS \n");
  //
  timer.start(t0f);
  ierr = l1.bfs_mesh_search(c0, querry);
  timer.stop();
  if(ierr != Point_locator::LOCATOR_SUCCESS) printf("fast failed ...\n");
  //
  if(perform_exact)
    {
      timer.start(t0e);
      ierr = l1.bruteforce_mesh_search(c0, querry);
      timer.stop();
      if(ierr != Point_locator::LOCATOR_SUCCESS) printf("exact failed ...\n");
    }
  //

  const int t1f = timer.new_activity("SECOND TO FIRST -- FAST");
  const int t1e = timer.new_activity("SECOND TO FIRST -- EXACT");
  fprintf(stderr, "*** SECOND TO FIRST PROJECTION STATUS \n");
  timer.start(t1f);
  ierr = l0.bfs_mesh_search(c1, querry);
  timer.stop();
  if(ierr != Point_locator::LOCATOR_SUCCESS) printf("fast failed ...\n");
  //
  if(perform_exact)
    {
      timer.start(t1e);
      ierr = l0.bruteforce_mesh_search(c1, querry);
      timer.stop();
      if(ierr != Point_locator::LOCATOR_SUCCESS) printf("exact failed ...\n");
    }
  //

  return 0;

} // End of main
