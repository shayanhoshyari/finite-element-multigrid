#ifndef BRIDSON_PCG_SOLVER_H
#define BRIDSON_PCG_SOLVER_H

// Implements PCG with Modified Incomplete Cholesky (0) preconditioner.
// PCGSolver<T> is the main class for setting up and solving a linear system.
// Note that this only handles symmetric positive (semi-)definite matrices,
// with guarantees made only for M-matrices (where off-diagonal entries are all
// non-positive, and row sums are non-negative).

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>

#include "blas_wrappers.h"
#include "icc.h"
#include "sparse_matrix.h"

namespace bridson
{


// ================================================================
//                     ABSTRACT PRECONDITIONER
// ================================================================
struct Preconditioner
{
  // Virtual destructor is very important. Otherwise, memory leak.
  virtual ~Preconditioner() {}
  // This is used so that we can make polymorphic copies
  virtual Preconditioner * clone() const = 0;
  //
  virtual void form(const FixedSparseMatrix & matrix) = 0;
  virtual void apply(const std::vector<double> & x, std::vector<double> & result) = 0;
};

// ================================================================
//                     NO PRECONDITIONER
// ================================================================
class Dummy_preconditioner : public Preconditioner
{
public:
  Preconditioner * clone() const { return new Dummy_preconditioner; }
  void form(const FixedSparseMatrix & matrix) {}
  void apply(const std::vector<double> & x, std::vector<double> & result) { result = x; }
};

// ================================================================
//                     JACOBI PRECONDITIONER
// ================================================================

class Jacobi_preconditioner : public Preconditioner
{
public:
  Preconditioner * clone() const { return new Jacobi_preconditioner; }
  void form(const FixedSparseMatrix & matrix)
  {
    _invdiag.resize(matrix.n);
    for(int i = 0; i < matrix.n; ++i)
      {
        _invdiag[i] = 1 / matrix.get_value(i, i);
      }
  }

  void apply(const std::vector<double> & x, std::vector<double> & result)
  {
    for(int i = 0; i < (int)_invdiag.size(); ++i)
      {
        result[i] = x[i] * _invdiag[i];
      }
  }

private:
  std::vector<double> _invdiag;
};

// ================================================================
//                     ICC0 PRECONDITIONER
// ================================================================
class ICC_preconditioner : public Preconditioner
{
public:
  Preconditioner * clone() const { return new ICC_preconditioner; }
  void form(const FixedSparseMatrix & matrix)
  { //
    factor_modified_incomplete_cholesky0(matrix, ic_factor);
  }
  void apply(const std::vector<double> & x, std::vector<double> & result)
  {
    solve_lower(ic_factor, x, result);
    solve_lower_transpose_in_place(ic_factor, result);
  }

private:
  SparseColumnLowerFactor ic_factor;
};


// ================================================================
//            SYMMETRIC GAUSS SEIDEL PRECONDITIONER
// ================================================================

//
// This is the symmetric Gauss-Sidel preconditioner.
// Non-symmetric versions will break the CG solver.
// Each iteration solves result = (D+R)^-1 D (D+L)^-1  x,
// in two steps:
// 1) y = (D+L)^-1  x
// 2) result = (D+R)^-1 D y
// https://math.stackexchange.com/questions/1484962/symmetric-gauss-seidel-iteration
//
// ASSUMPTION: columns are sorted for each row in the FixedSparseMatrix.
//
class Gauss_seidel_preconditioner : public Preconditioner
{
public:
  Preconditioner * clone() const { return new Gauss_seidel_preconditioner; }

  void form(const FixedSparseMatrix & matrix)
  {
    _invdiag.resize(matrix.n);
    _y.resize(matrix.n);
    _matrix = &matrix;
    for(int i = 0; i < matrix.n; ++i)
      {
        _invdiag[i] = 1 / matrix.get_value(i, i);
      }
  }

  void apply(const std::vector<double> & x, std::vector<double> & result)
  {
    double sum;

    for(int i = 0; i < (int)_invdiag.size(); ++i)
      {
        sum = 0;
        for(int flatid = _matrix->xadj[i]; flatid < _matrix->xadj[i + 1]; ++flatid)
          {
            const int j = _matrix->adj[flatid];
            if(j == i) break;
            assert(j < i);
            sum += _matrix->values[flatid] * _y[j];
          }
        _y[i] = (x[i] - sum) * _invdiag[i];
      }
    //
    for(int i = (int)_invdiag.size() - 1; i >= 0; --i)
      {
        sum = 0;
        for(int flatid = _matrix->xadj[i + 1] - 1; flatid >= _matrix->xadj[i]; --flatid)
          {
            const int j = _matrix->adj[flatid];
            if(j == i) break;
            assert(j > i);
            sum += _matrix->values[flatid] * result[j];
          }
        result[i] = _y[i] - sum * _invdiag[i];
      }
  }

private:
  std::vector<double> _invdiag, _y;
  const FixedSparseMatrix * _matrix;
};


// ================================================================
//       PRECONDITIONED CONJUGATE GRADIENT LINEAR SOLVER
// ================================================================

// Encapsulates the Conjugate Gradient algorithm.

struct PCGSolver
{
  PCGSolver()
  : preconditioner_memory_manager(new ICC_preconditioner)
  , preconditioner(preconditioner_memory_manager.get())
  {
    set_solver_parameters(1e-5, 100, 0.97, 0.25);
  }

  PCGSolver(Preconditioner * p, bool has_ownership = false)
  : preconditioner(p)
  {
    set_solver_parameters(1e-5, 100, 0.97, 0.25);
    if(has_ownership)
      {
        preconditioner_memory_manager.reset(preconditioner);
      }
  }

  void set_solver_parameters(double tolerance_factor_,
                             int max_iterations_,
                             double modified_incomplete_cholesky_parameter_ = 0.97,
                             double min_diagonal_ratio_ = 0.25)
  {
    tolerance_factor = tolerance_factor_;
    if(tolerance_factor < 1e-30) tolerance_factor = 1e-30;
    max_iterations = max_iterations_;
    modified_incomplete_cholesky_parameter = modified_incomplete_cholesky_parameter_;
    min_diagonal_ratio = min_diagonal_ratio_;
  }

  void set_max_iter(const int max_iterations_) { max_iterations = max_iterations_; }

  bool solve(const FixedSparseMatrix & matrix,
             const std::vector<double> & rhs,
             std::vector<double> & result,
             double & residual_out,
             int & iterations_out)
  {
    history.resize(0);
    history.reserve(max_iterations);
    unsigned int n = matrix.n;
    if(m.size() != n)
      {
        m.resize(n);
        s.resize(n);
        z.resize(n);
        r.resize(n);
      }
    zero(result);
    r = rhs;
    residual_out = BLAS::abs_max(r);
    history.push_back(residual_out);
    if(residual_out == 0)
      {
        iterations_out = 0;
        return true;
      }
    double tol = tolerance_factor * residual_out;

    preconditioner->form(matrix);
    preconditioner->apply(r, z);
    double rho = BLAS::dot(z, r);
    if(rho == 0 || rho != rho)
      {
        iterations_out = 0;
        return false;
      }

    s = z;
    int iteration;
    for(iteration = 0; iteration < max_iterations; ++iteration)
      {
        multiply(matrix, s, z);
        double alpha = rho / BLAS::dot(s, z);
        BLAS::add_scaled(alpha, s, result);
        BLAS::add_scaled(-alpha, z, r);
        residual_out = BLAS::abs_max(r);
        history.push_back(residual_out);
        if(residual_out <= tol)
          {
            iterations_out = iteration + 1;
            return true;
          }
        preconditioner->apply(r, z);
        double rho_new = BLAS::dot(z, r);
        double beta = rho_new / rho;
        BLAS::add_scaled(beta, s, z);
        s.swap(z); // s=beta*s+z
        rho = rho_new;
      }
    iterations_out = iteration;
    return false;
  }

public:
  std::vector<double> history;

protected:
  // internal structures
  std::unique_ptr<Preconditioner> preconditioner_memory_manager; // modified incomplete cholesky factor
  Preconditioner * preconditioner;
  std::vector<double> m, z, s, r; // temporary vectors for PCG
  // FixedSparseMatrix         fixed_matrix;                               // used within loop

  // parameters
  double tolerance_factor;
  int max_iterations;
  double modified_incomplete_cholesky_parameter;
  double min_diagonal_ratio;
};

} // end of namespace bridson

#endif
