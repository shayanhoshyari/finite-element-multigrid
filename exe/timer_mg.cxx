#include <cmath>
#include <fstream>

#include "../src/grid_io.hxx"
#include "../src/mg_tools.hxx"
#include "../src/poisson_equation.hxx"

#include "../src/bridson/pcg_solver.h"

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
}

int
main(int argc, char * argv[])
{
  FILE * fl = stderr;
  int use_pcg = 0;
  int max_mesh_index = 4;

  {
    for(int i = 0; i < argc; ++i)
      {
        auto option_present = [i, argc, argv](std::string name) -> bool {
          return (std::string(argv[i]) == name) && ((i + 1) < argc);
        };
        if(option_present("--maxindex")) sscanf(argv[i + 1], "%d", &max_mesh_index);
        if(option_present("--usepcg")) sscanf(argv[i + 1], "%d", &use_pcg);
      }
  }

  // =================================================
  mg::Context mg_context;
  Numeric_array rhs, solution;

  // =================================================
  Problem_data problem_data;
  Poisson_assembler assembler(problem_data);

  // =================================================
  {
    mg::Parameters mg_params;
    mg_params.cycle_type = mg::Parameters::CYCLE_TYPE_V;
    //
    mg_params.n_per_level_smoothing.resize(max_mesh_index);
    std::fill(mg_params.n_per_level_smoothing.begin(), mg_params.n_per_level_smoothing.end(), 3);
    mg_params.n_per_level_smoothing.back() = 500;
    //
    mg_context.set_parameters(mg_params);
  }

  // =================================================
  {
    std::vector<Grid_connectivity_context> grids(max_mesh_index);
    char mesh_file_name[1024];

    for(int i = 0; i < max_mesh_index; ++i)
      {
        sprintf(mesh_file_name, "meshes/squaremg.%d", max_mesh_index - i);
        Grid_IO(grids[i]).read_triangle(mesh_file_name);
        printf("Mesh level %d is %s \n", i, mesh_file_name);
      }

    mg_context.set_grids(grids);
  }

  // =================================================
  mg_context.set_up_discrete_equations(assembler);
  mg_context.set_up_transfer_info();
  mg_context.set_smoothers(bridson::Gauss_seidel_preconditioner());
  rhs.init(mg_context.ff[0][mg::ZERO].size());
  rhs.set_equal(mg_context.ff[0][mg::ZERO]);
  solution.init(mg_context.ff[0][mg::ZERO].size());

  // ================================================
  std::vector<double> history;
  int n_cg_iterations;
  double cg_residual;
  if(!use_pcg)
    {
      mg_context.set_solver_parameters(1e-15, 15);
      mg_context.solve(rhs, solution, cg_residual, n_cg_iterations);
      history = mg_context.history;
    }
  else
    {
      // In this case the RHS has to be provided separately
      bridson::PCGSolver pcg_solver(&mg_context, false);
      //
      // bridson::PCGSolver pcg_solver(new bridson::Gauss_seidel_preconditioner, true);
      //
      pcg_solver.set_solver_parameters(1e-15, 10000);

      pcg_solver.solve(mg_context.AA[0], rhs, solution, cg_residual, n_cg_iterations);
      history = pcg_solver.history;
    }
  fprintf(fl, "solver: residual: %g n_iter: %d \n", cg_residual, n_cg_iterations);
  for(unsigned i = 0; i < history.size(); ++i)
    {
      printf("Iteration %5d, residual: %15.8g \n", (int)i, history[i]);
    }

  // ================================================

  // ================================================
  const int n_dofs = assembler.get_n_dof(mg_context.grids[0]);
  Numeric_array dof_exact(n_dofs);
  Numeric_array dof_diff(n_dofs);
  Numeric_array & dof_numeric = solution;

  assembler.find_exact_dofs(mg_context.grids[0], dof_exact);
  dof_diff.set_equal(dof_exact);
  dof_diff.add_scaled(-1, dof_numeric);
  fprintf(fl,
          "Solution errors: %20.10g %20.10g %20.10g\n",
          dof_diff.norm_l1() / n_dofs,
          dof_diff.norm_l2() / sqrt(n_dofs),
          dof_diff.norm_linf());

  {
    // For octave actually. The whole matlab gui is so
    // sluggy for this quick checks.
    // std::ofstream os("formatlab.m");
    // lhs.write_matlab(os, "lhs");
    // rhs.write_matlab(os, "rhs");
    // dof_numeric.write_matlab(os, "dof_numeric");
    // dof_exact.write_matlab(os, "dof_exact");
  }

  Grid_IO io(mg_context.grids[0]);
  FILE * vtk_fl = fopen("__poisson_test.2.vtk", "w");
  io.write_vtk(vtk_fl);
  io.write_vtk_vert_header(vtk_fl);
  io.write_vtk_data(vtk_fl, dof_exact, "u_exact");
  io.write_vtk_data(vtk_fl, dof_numeric, "u_numeric");
  io.write_vtk_data(vtk_fl, dof_diff, "error");
  fclose(vtk_fl);


  return 0;
}
