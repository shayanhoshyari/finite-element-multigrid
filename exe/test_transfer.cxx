#include <cmath>
#include <fstream>

#include "../src/grid_io.hxx"
#include "../src/grid_iterators.hxx"
#include "../src/poisson_equation.hxx"

#include "../src/mg_tools.hxx"

namespace tests
{

double
func(double x, double y)
{
  constexpr double pi = 3.141592653589793238462643383279502884;
  return sin(pi * x) * sin(2 * pi * y);
}

double
test_transfer_once(FILE * fl, std::string mcoarse, std::string mfine)
{
  // READ MESHES
  Grid_connectivity_context gfine, gcoarse;
  Grid_iterators giterfine(gfine), gitercoarse(gcoarse);
  Grid_IO(gfine).read_triangle(mfine);
  Grid_IO(gcoarse).read_triangle(mcoarse);

  // Create the DOF vectors
  Numeric_array ffine(gfine.n_vertices());
  Numeric_array fcoarse(gcoarse.n_vertices());
  Numeric_array fc2f(gfine.n_vertices());
  Numeric_array error(gfine.n_vertices());
  //
  Numeric_array fc2f2c(gcoarse.n_vertices());
  Numeric_array errorc(gcoarse.n_vertices());

  // Find locator
  Point_locator::Grid_query_out transfer_data;
  const int ierr = Point_locator(gcoarse).bfs_mesh_search(gfine, transfer_data);
  assert(ierr == Point_locator::LOCATOR_SUCCESS);

  // Find function three ways
  // On fine
  // On fine transfered to coarse
  // On coarse
  for(int i = 0; i < giterfine.n_vertices(); ++i)
    {
      Grid_iterators::Vertex v = giterfine.vertex_iterator(i);
      const double f = func(v.x(), v.y());
      ffine[i] = f;
    }
  for(int i = 0; i < gitercoarse.n_vertices(); ++i)
    {
      Grid_iterators::Vertex v = gitercoarse.vertex_iterator(i);
      const double f = func(v.x(), v.y());
      fcoarse[i] = f;
    }
  mg::Context::interpolate_impl(transfer_data, gcoarse, fcoarse, gfine, fc2f);
  error.set_equal(ffine);
  error.add_scaled(-1, fc2f);

  mg::Context::interpolate_transpose_impl(transfer_data, gfine, fc2f, gcoarse, fc2f2c);
  errorc.set_equal(fcoarse);
  errorc.add_scaled(-1, fc2f2c);

  fprintf(fl, "%8s to %8s: %10.8g ", mcoarse.c_str(), mfine.c_str(), error.norm_linf());

  {
    Grid_IO io(gfine);
    FILE * vtk_fl = fopen("_mesh_fine.vtk", "w");
    io.write_vtk(vtk_fl);
    io.write_vtk_vert_header(vtk_fl);
    io.write_vtk_data(vtk_fl, ffine, "u_exact");
    io.write_vtk_data(vtk_fl, fc2f,  "u_interpolate");
    io.write_vtk_data(vtk_fl, error, "error");
    fclose(vtk_fl);
  }
  {
    Grid_IO io(gcoarse);
    FILE * vtk_fl = fopen("_mesh_coarse.vtk", "w");
    io.write_vtk(vtk_fl);
    io.write_vtk_vert_header(vtk_fl);
    io.write_vtk_data(vtk_fl, fcoarse, "u_exact");
    io.write_vtk_data(vtk_fl, fc2f2c, "u_transfer");
    io.write_vtk_data(vtk_fl, errorc, "error");
    fclose(vtk_fl);
  }

  return error.norm_linf();
}

void
test_transfer(FILE * fl)
{
  const double e1 = test_transfer_once(fl, "meshes/squaremg.1", "meshes/squaremg.2");
  fprintf(fl,"\n");
  const double e2 = test_transfer_once(fl, "meshes/squaremg.2", "meshes/squaremg.3");
  fprintf(fl," ratio: %g \n", e1/e2);
  const double e3 = test_transfer_once(fl, "meshes/squaremg.3", "meshes/squaremg.4");
  fprintf(fl," ratio: %g \n", e2/e3);
}

} // end of tests