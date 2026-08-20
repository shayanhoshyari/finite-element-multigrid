
#include "../src/bridson/pcg_solver.h"
#include "../src/fake_fprintf.hxx"

#include "../src/numeric_array.hxx"

#include <cassert>

namespace tests
{
bool
numerics(FILE * fl)
{
  using namespace bridson;
  File_stream flstream(fl);

  Numeric_array x(4);
  x.set_equal(1);
  fprintf(fl, "%s", x.to_string().c_str());
  x.add_scaled(1, x);
  fprintf(fl, "%s", x.to_string().c_str());


  // Create a sparse system
  std::vector<int> adj, xadj;
  adj = {0, 2, 3, 1, 2, 0, 1, 2, 0, 3};
  xadj = {0, 3, 5, 8, 10};
  int n_iterations;
  double residual;

  FixedSparseMatrix sp;
  Numeric_array rhs(4), soln(4), sole(4);
  rhs.set_equal({40, 32, 27, 86});
  sole.set_equal({4, 3, 2, 7});


  {
    // Set the lhs directly
    sp.set_adjacency(xadj, adj);
    sp.values = {2, 2, 4, 10, 1, 2, 1, 8, 4, 10};

    PCGSolver cg(new bridson::Dummy_preconditioner, true);
    cg.set_solver_parameters(1e-30, 1000);
    cg.solve(sp, rhs, soln, residual, n_iterations);

    flstream << "======== Setting vals directly\n";
    fprintf(fl, "n_iterations: %d, residual: %lf\n", n_iterations, residual);
    sp.write_matlab(flstream, "lhs");
    soln.write_matlab(flstream, "soln");
    sole.write_matlab(flstream, "sole");
    rhs.write_matlab(flstream, "rhs");
  }

  {
    // Use set values
    sp.set_zero();
    sp.set_value(0, 0, 2);
    sp.set_value(0, 2, 2);
    sp.set_value(0, 3, 4);
    //
    sp.set_value(1, 1, 10);
    sp.set_value(1, 2, 1);
    //
    sp.set_value(2, 0, 2);
    sp.set_value(2, 1, 1);
    sp.set_value(2, 2, 8);
    //
    sp.set_value(3, 0, 4);
    sp.set_value(3, 3, 10);

    PCGSolver cg(new bridson::ICC_preconditioner, true);
    cg.set_solver_parameters(1e-30, 1000);
    cg.solve(sp, rhs, soln, residual, n_iterations);

    flstream << "======== Setting vals using set_value\n";
    fprintf(fl, "n_iterations: %d, residual: %lf\n", n_iterations, residual);
    sp.write_matlab(flstream, "lhs");
    soln.write_matlab(flstream, "soln");
    sole.write_matlab(flstream, "sole");
    rhs.write_matlab(flstream, "rhs");
  }

  {
    // Use add values
    sp.set_zero();
    sp.add_values(Vec2i(0,2), Vec2i(0,2), Vec4d(2,2,2,8) );
    sp.add_values(Vec2i(2,1), Vec1i(1), Vec2d(1,10) );
    sp.add_values(Vec2i(0,3), Vec1i(3), Vec2d(4,10) );
    sp.add_values(Vec1i(1), Vec1i(2), Vec1d(1) );
    sp.add_values(Vec1i(3), Vec1i(0), Vec1d(4) );

    PCGSolver cg(new bridson::Gauss_seidel_preconditioner, true);
    cg.set_solver_parameters(1e-30, 1000);
    cg.solve(sp, rhs, soln, residual, n_iterations);

    flstream << "======== Setting vals using set_value\n";
    fprintf(fl, "n_iterations: %d, residual: %lf\n", n_iterations, residual);
    sp.write_matlab(flstream, "lhs");
    soln.write_matlab(flstream, "soln");
    sole.write_matlab(flstream, "sole");
    rhs.write_matlab(flstream, "rhs");
  }


  return true;
}
}