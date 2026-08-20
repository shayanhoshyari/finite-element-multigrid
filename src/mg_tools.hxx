#if !defined(MG_TOOLS_IS_INCLUDED)
#define MG_TOOLS_IS_INCLUDED

#include <array>
#include <memory>
#include <vector>

#include "bridson/pcg_solver.h"
#include "bridson/sparse_matrix.h"
#include "grid_connectivity_context.hxx"
#include "numeric_array.hxx"
#include "point_locator.hxx"
#include "poisson_equation.hxx"

namespace mg
{

// Each vector has two indices, e.g., uu[level_id][buffer_id].
// Use numbers for the buffer index. Otherwise, the code becomes
// quite confusing.
static constexpr int n_buffers = 1;
enum Numbers
{
  ZERO = 0,
  ONE,
  TWO,
  THREE,
  FOUR,
  FIVE,
  SIX,
  SEVEN,
  EIGHT,
  NINE,
  TEN
};

struct Parameters
{
  enum Cycle_type
  {
    CYCLE_TYPE_V = 0,
    CYCLE_TYPE_W,
  };

  Cycle_type cycle_type;
  double relaxation_factor;
  std::vector<int> n_per_level_smoothing; // post and pre

  int n_levels() const { return (int)n_per_level_smoothing.size(); }
};

//
// I visualize the multigrid hierarchy as follows:
//
//         Level n_levels-1    Coarsest
//           /
//         ...
//         /
//      Level 3
//       /
//    Level 1
//     /
// Level 0                     Finest (on which we solve for the solution)
//
//  transfer_info[i] (size n_levels-2): location of level i on i+1
//  grids,
class Context : public ::bridson::Preconditioner
{
public:
  // ==================================================
  // Public members
  // ==================================================
  typedef Point_locator::Grid_query_out Transfer_info;
  typedef bridson::FixedSparseMatrix Sparse_matrix;
  typedef bridson::Preconditioner Preconditioner;

  // Parameters
  Parameters parameters;

  // Grids and transfer
  std::vector<Grid_connectivity_context> grids;
  std::vector<Transfer_info> transfer_infos;

  // Numerics
  std::vector<Sparse_matrix> AA;                          // left hand side
  std::vector<std::unique_ptr<Preconditioner>> smoothers; // smoothers
  std::vector<std::array<Numeric_array, n_buffers>> uu;   // solution
  std::vector<std::array<Numeric_array, n_buffers>> rr;   // residual
  std::vector<std::array<Numeric_array, n_buffers>> du;   // solution update
  std::vector<std::array<Numeric_array, n_buffers>> ff;   // right hand side

  // History of residual
  std::vector<double> history;

  // ==================================================
  // Constructor
  // ==================================================
  Context();

  // ==================================================
  // Basic info
  // ==================================================

  int n_levels() const
  {
    assert(_is_set_parameters_called);
    return parameters.n_levels();
  }
  int n_smoothing(const int level) const
  {
    assert(_is_set_parameters_called);
    assert(level < n_levels());
    return parameters.n_per_level_smoothing[level];
  }
  void verify() const;

  // ==================================================
  // Preparing the system
  // ==================================================

  void set_parameters(const Parameters &);
  void set_grids(std::vector<Grid_connectivity_context> & grids_in); // Swaps the grids
  void set_up_transfer_info();
  void set_smoothers(std::vector<Preconditioner *> & precons_in); // takes ownership
  void set_smoothers(const Preconditioner & precon_in);           // only clones
  void set_up_discrete_equations(const Poisson_assembler &);

  // ==================================================
  // Main multigrid actions
  // ==================================================

  // Prolongation would be  interpolation
  void prolongate_dofs(const int level_from, const int level_to, const Numeric_array & in, Numeric_array & out);
  // Restriction would be the transpose of interpolation
  void restrict_dofs(const int level_from, const int level_to, const Numeric_array & in, Numeric_array & out);
  // For internal calls
  // TODO: use a different impl for Transfer_info and make this faster.
  static void interpolate_impl(const Transfer_info & tr,
                               const Grid_connectivity_context & gin,
                               const Numeric_array & in,
                               const Grid_connectivity_context & gout,
                               Numeric_array & out);
  static void interpolate_transpose_impl(const Transfer_info & tr,
                                         const Grid_connectivity_context & gin,
                                         const Numeric_array & in,
                                         const Grid_connectivity_context & gout,
                                         Numeric_array & out);

  // Add eval_residual
  static double eval_residual(const Sparse_matrix & A,
                              const Numeric_array & u,
                              const Numeric_array & f,
                              Numeric_array & r);

  // Assumes we have the system A[i]
  // returns previous residual
  double smooth(const int level, const int buffer_id, const int n_times);

  // Cycle
  // returns previous residual
  double cycle(const int level, const int buffer_id);

  // ==================================================
  // Preconditioner Interface
  // ==================================================
  bridson::Preconditioner * clone() const
  {
    assert(0 && "Not implemented");
    return NULL;
  }
  //
  void form(const bridson::FixedSparseMatrix & matrix) {}
  void apply(const std::vector<double> & x, std::vector<double> & result)
  {
    ff[0][ZERO].set_equal(x);
    uu[0][ZERO].set_zero();
    this->cycle(0, ZERO);
    std::copy(uu[0][ZERO].begin(), uu[0][ZERO].end(), result.begin());
  }

  // ==================================================
  // Independent solver interface
  // ==================================================
  bool solve(const std::vector<double> & rhs, std::vector<double> & result, double & residual_out, int & iterations_out)
  {
    bool success = false;

    ff[0][ZERO].set_equal(rhs);
    uu[0][ZERO].set_zero();
    history.resize(0);
    history.reserve(_max_iterations);
    for(int i = 0; i < _max_iterations; ++i)
      {
        residual_out = this->cycle(0, ZERO);
        history.push_back( residual_out );
        if( residual_out < _tolerance_factor)
        {
          success = true;
          break;
        }
      }

    iterations_out = (int)history.size();
    std::copy(uu[0][ZERO].begin(), uu[0][ZERO].end(), result.begin());

    return success;
  }

  void set_solver_parameters(double tolerance_factor_,int max_iterations_)
  {
    _max_iterations = max_iterations_;
    _tolerance_factor = tolerance_factor_;
  }

private:
  bool _is_set_smoothers_called;
  bool _is_set_grids_called;
  bool _is_set_parameters_called;
  bool _is_set_up_transfer_info_called;
  bool _is_set_up_discrete_equations_called;

  // solver parameters
  double _tolerance_factor;
  int _max_iterations;

  // Don't copy this object around
  Context & operator=(const Context &) = delete;
  Context(const Context &) = delete;
};



} // End of namespace mg

#endif /* MG_TOOLS_IS_INCLUDED */