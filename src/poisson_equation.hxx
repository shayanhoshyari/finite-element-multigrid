#ifndef POISSON_EQUATION_IS_INCLUDED
#define POISSON_EQUATION_IS_INCLUDED

#include "bridson/sparse_matrix.h"
#include "bridson/vec.h"
#include "numeric_array.hxx"
#include "grid_connectivity_context.hxx"

/****************************************************************************
 *                        Problem data                                      *
 ****************************************************************************/
// -div( K grad(u)) = f
class Poisson_problem_data
{
public:
  virtual ~Poisson_problem_data() {}

  //
  // Exact solution if available.
  //
  virtual bool has_exact_solution() const { return false; }
  virtual double exact_solution(const bridson::Vec2d &) const
  {
    assert(0);
    return 0;
  }

  //
  // Exact derivative if available.
  //
  virtual bool has_exact_derivative() const { return false; }
  virtual bridson::Vec4d exact_derivative(const bridson::Vec2d &) const
  {
    assert(0);
    return bridson::Vec4d();
  }

  //
  // The K tensor as in -div( K grad(phi) ) = f
  // The stiffness tensor is 2*2 for poisson
  // In general it is a fourth order tensor n_unknown*n_dim*n_unknown*n_dim or it can be
  // flattened into a (n_unknown*n_dim)*(n_unknown*n_dim) matrix.
  //
  virtual bridson::Vec4d stiffness() const = 0;

  //
  // Right hand side of the equations.
  //
  virtual double right_hand_side(const bridson::Vec2d &) const = 0;

  //
  // Boundary conditions
  // Only supported type is ROBIN, i.e., q(x)u + \nabla u \cdot n = t(x)
  //
  virtual void boundary_condition_terms(const bridson::Vec2d &, const int btag, double & q, double & t) const = 0;
};

/****************************************************************************
 *                             Assembler                                    *
 ****************************************************************************/

//
// Thoughts: in case we want to handle a wider range of basis functions or
// equations, rather than a poisson_assembler, the functionality should be decomposed:
//
// DOF handler
//    gets mesh + FE type -> spits non-zero structure
//    get element + FE type + dof -> says where it is non zero
//    etc.
//
// Assembling: Allow user to specify via their own version of XXX_assembler
//

// TODO: test quadrature
//
//
class Poisson_assembler
{
public:
  Poisson_assembler(const Poisson_problem_data & problem_data_in)
  : _problem_data(&problem_data_in)
  {
  }

  /*
   * Number of DOF for a mesh
   */
  int get_n_dof(const Grid_connectivity_context & grid) const;

  /**
   * Assemble the system matrix and right-hand side vector.
   * Follows the sparse row format in FixedSparseMatrix or the METIS graph partitioner
   * adj[ xadj[i]:xadj[i+1] ] are the neighbours of the DOF i.
   */
  void get_nonzero_structure(const Grid_connectivity_context & grid,
                             std::vector<int> & xadj,
                             std::vector<int> & adj) const;
  void get_nonzero_structure(const Grid_connectivity_context & grid, bridson::FixedSparseMatrix & lhs) const;

  /**
   * Assembles the system of equations
   * Does not
   */
  void assemble_system(const Grid_connectivity_context & grid,
                       bridson::FixedSparseMatrix & lhs,
                       Numeric_array & rhs) const;


  //
  // For matrix free smoother
  // later
  // void matrix_free_lhs_product(const Grid_connectivity_context & grid,
  //                             const Numeric_array & in,
  //                             Numeric_array & out) const;


  //
  // Project the exact solution on the dofs
  // Uses interpolation
  //
  void find_exact_dofs(const Grid_connectivity_context & grid, Numeric_array & dofs) const;

  const Poisson_problem_data & problem_data() const { return *_problem_data; }

private:
  const Poisson_problem_data * _problem_data;
};

#endif /* POISSON_EQUATION_IS_INCLUDED */
