#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "../src/bridson/pcg_solver.h"
#include "../src/grid_io.hxx"
#include "../src/mg_tools.hxx"
#include "../src/poisson_equation.hxx"

namespace
{
class Problem_data : public Poisson_problem_data
{
  static constexpr double pi = 3.141592653589793238462643383279502884;
  static constexpr double pi2 = pi * pi;

  bool has_exact_solution() const { return true; }
  double exact_solution(const bridson::Vec2d & p) const { return std::sin(pi * p(0)) * std::sin(2 * pi * p(1)); }
  bridson::Vec4d stiffness() const { return bridson::Vec4d(1, 0, 0, 1); }
  double right_hand_side(const bridson::Vec2d & p) const { return 5 * pi2 * exact_solution(p); }
  void boundary_condition_terms(const bridson::Vec2d & p, const int btag, double & q, double & t) const
  {
    const double beta = 1e10;
    const bool valid = btag == 101 || btag == 102 || btag == 103 || btag == 104;
    if(!valid)
      {
        std::fprintf(stderr, "Invalid boundary tag: %d\n", btag);
        std::abort();
      }
    q = beta;
    t = beta * exact_solution(p);
  }
};

struct Options
{
  std::string method = "mg-v";
  int fine_level = 4;
  int max_iterations = 100;
  int smoothing = 3;
  int coarse_smoothing = 500;
  double tolerance = 1e-12;
  bool print_history = false;
  std::string field_output;
};

void
usage(const char * name)
{
  std::fprintf(stderr,
               "Usage: %s [--method mg-v|mg-w|gs|pcg-none|pcg-gs|pcg-mg-v|pcg-mg-w]\n"
               "          [--fine N] [--max-iterations N] [--smoothing N]\n"
               "          [--coarse-smoothing N] [--tolerance X] [--history]\n"
               "          [--field-output FILE]\n",
               name);
}

Options
parse_options(int argc, char ** argv)
{
  Options options;
  for(int i = 1; i < argc; ++i)
    {
      const std::string argument(argv[i]);
      if(argument == "--history")
        options.print_history = true;
      else if(argument == "--help")
        {
          usage(argv[0]);
          std::exit(0);
        }
      else
        {
          if(i + 1 == argc)
            {
              usage(argv[0]);
              std::exit(2);
            }
          const std::string value(argv[++i]);
          if(argument == "--method")
            options.method = value;
          else if(argument == "--field-output")
            options.field_output = value;
          else if(argument == "--fine")
            options.fine_level = std::atoi(value.c_str());
          else if(argument == "--max-iterations")
            options.max_iterations = std::atoi(value.c_str());
          else if(argument == "--smoothing")
            options.smoothing = std::atoi(value.c_str());
          else if(argument == "--coarse-smoothing")
            options.coarse_smoothing = std::atoi(value.c_str());
          else if(argument == "--tolerance")
            options.tolerance = std::atof(value.c_str());
          else
            {
              usage(argv[0]);
              std::exit(2);
            }
        }
    }
  return options;
}

bool
is_multilevel(const std::string & method)
{
  return method == "mg-v" || method == "mg-w" || method == "pcg-mg-v" || method == "pcg-mg-w";
}

bool
is_w_cycle(const std::string & method)
{
  return method == "mg-w" || method == "pcg-mg-w";
}

bool
is_pcg(const std::string & method)
{
  return method.compare(0, 3, "pcg") == 0;
}
} // namespace

int
main(int argc, char ** argv)
{
  const Options options = parse_options(argc, argv);
  const std::vector<std::string> valid_methods = {"mg-v", "mg-w", "gs", "pcg-none", "pcg-gs", "pcg-mg-v", "pcg-mg-w"};
  if(std::find(valid_methods.begin(), valid_methods.end(), options.method) == valid_methods.end() ||
     options.fine_level < 1)
    {
      usage(argv[0]);
      return 2;
    }

  const int n_levels = is_multilevel(options.method) ? options.fine_level : 1;
  mg::Parameters parameters;
  parameters.cycle_type = is_w_cycle(options.method) ? mg::Parameters::CYCLE_TYPE_W : mg::Parameters::CYCLE_TYPE_V;
  parameters.relaxation_factor = 1;
  parameters.n_per_level_smoothing.assign(n_levels, options.smoothing);
  if(options.method == "gs")
    parameters.n_per_level_smoothing.back() = 1;
  else
    parameters.n_per_level_smoothing.back() = options.coarse_smoothing;

  mg::Context context;
  context.set_parameters(parameters);

  std::vector<Grid_connectivity_context> grids(n_levels);
  for(int level = 0; level < n_levels; ++level)
    {
      const int mesh_index = is_multilevel(options.method) ? options.fine_level - level : options.fine_level;
      char filename[256];
      std::snprintf(filename, sizeof(filename), "meshes/squaremg.%d", mesh_index);
      Grid_IO(grids[level]).read_triangle(filename);
    }

  const int vertices = grids[0].n_vertices();
  const int triangles = grids[0].n_real_triangles();
  context.set_grids(grids);

  Problem_data problem_data;
  Poisson_assembler assembler(problem_data);
  const auto setup_start = std::chrono::steady_clock::now();
  context.set_up_discrete_equations(assembler);
  if(n_levels > 1)
    context.set_up_transfer_info();
  else
    {
      // A zero-transfer single-level context still needs this lifecycle flag.
      context.set_up_transfer_info();
    }
  context.set_smoothers(bridson::Gauss_seidel_preconditioner());
  const auto setup_end = std::chrono::steady_clock::now();

  Numeric_array rhs(context.ff[0][mg::ZERO].size());
  Numeric_array solution(rhs.size());
  rhs.set_equal(context.ff[0][mg::ZERO]);

  bool success = false;
  double residual = 0;
  int iterations = 0;
  std::vector<double> history;
  const auto solve_start = std::chrono::steady_clock::now();
  if(!is_pcg(options.method))
    {
      context.set_solver_parameters(options.tolerance, options.max_iterations);
      success = context.solve(rhs, solution, residual, iterations);
      history = context.history;
    }
  else
    {
      bridson::Preconditioner * preconditioner = NULL;
      bool owns_preconditioner = true;
      if(options.method == "pcg-none")
        preconditioner = new bridson::Dummy_preconditioner;
      else if(options.method == "pcg-gs")
        preconditioner = new bridson::Gauss_seidel_preconditioner;
      else
        {
          preconditioner = &context;
          owns_preconditioner = false;
        }
      bridson::PCGSolver solver(preconditioner, owns_preconditioner);
      solver.set_solver_parameters(options.tolerance, options.max_iterations);
      success = solver.solve(context.AA[0], rhs, solution, residual, iterations);
      history = solver.history;
    }
  const auto solve_end = std::chrono::steady_clock::now();

  Numeric_array exact(solution.size());
  Numeric_array error(solution.size());
  assembler.find_exact_dofs(context.grids[0], exact);
  error.set_equal(exact);
  error.add_scaled(-1, solution);

  if(!options.field_output.empty())
    {
      Numeric_array final_residual(rhs.size());
      mg::Context::eval_residual(context.AA[0], solution, rhs, final_residual);
      FILE * field_file = std::fopen(options.field_output.c_str(), "w");
      if(!field_file)
        {
          std::fprintf(stderr, "Could not open field output: %s\n", options.field_output.c_str());
          return 2;
        }
      std::fprintf(field_file, "vertex,x,y,solution,residual\n");
      const std::vector<Grid_connectivity_context::Vertex> & vertices_data = context.grids[0].vertices_data();
      for(int vertex = 0; vertex < vertices; ++vertex)
        std::fprintf(field_file,
                     "%d,%.17g,%.17g,%.17g,%.17g\n",
                     vertex,
                     vertices_data[vertex].xy[0],
                     vertices_data[vertex].xy[1],
                     solution[vertex],
                     final_residual[vertex]);
      std::fclose(field_file);
    }

  const std::chrono::duration<double> setup_seconds = setup_end - setup_start;
  const std::chrono::duration<double> solve_seconds = solve_end - solve_start;
  std::printf("RESULT method=%s fine=%d levels=%d vertices=%d triangles=%d success=%d iterations=%d "
              "residual=%.17g setup_seconds=%.9g solve_seconds=%.9g l1=%.17g l2=%.17g linf=%.17g\n",
              options.method.c_str(),
              options.fine_level,
              n_levels,
              vertices,
              triangles,
              success ? 1 : 0,
              iterations,
              residual,
              setup_seconds.count(),
              solve_seconds.count(),
              error.norm_l1() / vertices,
              error.norm_l2() / std::sqrt(vertices),
              error.norm_linf());
  if(options.print_history)
    for(std::size_t i = 0; i < history.size(); ++i)
      std::printf("HISTORY method=%s fine=%d iteration=%zu residual=%.17g\n",
                  options.method.c_str(),
                  options.fine_level,
                  i,
                  history[i]);

  return success ? 0 : 1;
}
