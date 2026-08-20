#include <cmath>
#include <fstream>

#include "../src/grid_io.hxx"
#include "../src/poisson_equation.hxx"

#include "../src/bridson/pcg_solver.h"

namespace tests
{

// Anonymus namespace to prevent name clash for Problem_data.
namespace
{

class Problem_data : public Poisson_problem_data
{
  static constexpr double pi = 3.141592653589793238462643383279502884;
  static constexpr double pi2 = pi * pi;

  bool has_exact_solution() const { return true; }
  double exact_solution(const bridson::Vec2d & p) const
  {

    const double x = p(0);
    const double y = p(1);
    return sin(pi * x) * sin(2 * pi * y);
  }

  bridson::Vec4d stiffness() const { return bridson::Vec4d(1, 0, 0, 1); }

  double right_hand_side(const bridson::Vec2d & p) const
  {
    const double sole = exact_solution(p);
    const double d2udx2 = -pi2 * sole;
    const double d2udy2 = -4 * pi2 * sole;
    return -(d2udx2 + d2udy2);
  }

  void boundary_condition_terms(const bridson::Vec2d & p, const int btag, double & q, double & t) const
  {
    const double beta = 1e10;
    const bool is_valid_tag = ((btag == 101) || (btag == 102) || (btag == 103) || (btag == 104));
    if(!is_valid_tag)
      {
        printf("Wrong tag %d \n", btag);
        assert(0 && "Wrong tag");
        throw;
      }

    q = beta;
    t = beta * exact_solution(p);
  }
}; // End of problem data

} // End of anonymus namespace

// SQUARE 4 (963 triangles)
// PCG solver: residual: 1.11081e-31 n_iter: 91
// Solution errors: 0.001186 0.001706 0.006713
//
// SQUARE 5 (15500 triangles -> h / 4 -> error / 16)
// PCG solver: residual: 8.73598e-33 n_iter: 358
// Solution errors: 0.000073 0.000102 0.000715
void
test_poisson(FILE * fl)
{
  Grid_connectivity_context g;
  Grid_IO io(g);

  io.read_triangle("meshes/square.4");
  io.write_vtk("__poisson_test.1.vtk");

  Problem_data problem_data;
  Poisson_assembler assembler(problem_data);

  const int n_dofs = assembler.get_n_dof(g);
  //
  bridson::FixedSparseMatrix lhs;
  assembler.get_nonzero_structure(g, lhs);
  //
  Numeric_array rhs(n_dofs);
  Numeric_array dof_numeric(n_dofs);
  Numeric_array dof_exact(n_dofs);
  Numeric_array dof_diff(n_dofs);

  assembler.assemble_system(g, lhs, rhs);

  bridson::PCGSolver pcg_solver(new bridson::Gauss_seidel_preconditioner, true);
  pcg_solver.set_solver_parameters(1e-30, 10000);

  int n_cg_iterations;
  double cg_residual;
  pcg_solver.solve(lhs, rhs, dof_numeric, cg_residual, n_cg_iterations);
  fprintf(fl, "PCG solver: residual: %g n_iter: %d \n", cg_residual, n_cg_iterations);

  assembler.find_exact_dofs(g, dof_exact);
  dof_diff.set_equal(dof_exact);
  dof_diff.add_scaled(-1, dof_numeric);
  fprintf(fl,
          "Solution errors: %lf %lf %lf\n",
          dof_diff.norm_l1() / n_dofs,
          dof_diff.norm_l2() / sqrt(n_dofs),
          dof_diff.norm_linf());

  {
    // For octave actually. The whole matlab gui is so
    // sluggy for this quick checks.
    std::ofstream os("formatlab.m");
    lhs.write_matlab(os, "lhs");
    rhs.write_matlab(os, "rhs");
    dof_numeric.write_matlab(os, "dof_numeric");
    dof_exact.write_matlab(os, "dof_exact");
  }

  FILE * vtk_fl = fopen("__poisson_test.2.vtk", "w");
  io.write_vtk(vtk_fl);
  io.write_vtk_vert_header(vtk_fl);
  io.write_vtk_data(vtk_fl, dof_exact, "u_exact");
  io.write_vtk_data(vtk_fl, dof_numeric, "u_numeric");
  io.write_vtk_data(vtk_fl, dof_diff, "error");
  fclose(vtk_fl);
}
}