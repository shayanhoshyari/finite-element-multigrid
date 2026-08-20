#ifndef BRIDSON_ICC_H
#define BRIDSON_ICC_H

// Implements PCG with Modified Incomplete Cholesky (0) preconditioner.
// PCGSolver<T> is the main class for setting up and solving a linear system.
// Note that this only handles symmetric positive (semi-)definite matrices,
// with guarantees made only for M-matrices (where off-diagonal entries are all
// non-positive, and row sums are non-negative).

#include "blas_wrappers.h"
#include "sparse_matrix.h"
#include <cmath>

namespace bridson
{

//============================================================================
// A simple compressed sparse column data structure (with separate diagonal)
// for lower triangular matrices

struct SparseColumnLowerFactor
{
  unsigned int n;
  std::vector<double> invdiag;        // reciprocals of diagonal elements
  std::vector<double> value;          // values below the diagonal, listed column by column
  std::vector<unsigned int> rowindex; // a list of all row indices, for each column in turn
  std::vector<unsigned int>
    colstart;                // where each column begins in rowindex (plus an extra entry at the end, of #nonzeros)
  std::vector<double> adiag; // just used in factorization: minimum "safe" diagonal entry allowed

  explicit SparseColumnLowerFactor(unsigned int n_ = 0)
  : n(n_)
  , invdiag(n_)
  , colstart(n_ + 1)
  , adiag(n_)
  {
  }

  void clear(void)
  {
    n = 0;
    invdiag.clear();
    value.clear();
    rowindex.clear();
    colstart.clear();
    adiag.clear();
  }

  void resize(unsigned int n_)
  {
    n = n_;
    invdiag.resize(n);
    colstart.resize(n + 1);
    adiag.resize(n);
  }

  void write_matlab(std::ostream & output, const char * variable_name)
  {
    output << variable_name << "=sparse([";
    for(unsigned int i = 0; i < n; ++i)
      {
        output << " " << i + 1;
        for(unsigned int j = colstart[i]; j < colstart[i + 1]; ++j)
          {
            output << " " << rowindex[j] + 1;
          }
      }
    output << "],...\n  [";
    for(unsigned int i = 0; i < n; ++i)
      {
        output << " " << i + 1;
        for(unsigned int j = colstart[i]; j < colstart[i + 1]; ++j)
          {
            output << " " << i + 1;
          }
      }
    output << "],...\n  [";
    for(unsigned int i = 0; i < n; ++i)
      {
        output << " " << (invdiag[i] != 0 ? 1 / invdiag[i] : 0);
        for(unsigned int j = colstart[i]; j < colstart[i + 1]; ++j)
          {
            output << " " << value[j];
          }
      }
    output << "], " << n << ", " << n << ");" << std::endl;
  }
};

//============================================================================
// Incomplete Cholesky factorization, level zero, with option for modified version.
// Set modification_parameter between zero (regular incomplete Cholesky) and
// one (fully modified version), with values close to one usually giving the best
// results. The min_diagonal_ratio parameter is used to detect and correct
// problems in factorization: if a pivot is this much less than the diagonal
// entry from the original matrix, the original matrix entry is used instead.

inline void
factor_modified_incomplete_cholesky0(const FixedSparseMatrix & matrix,
                                     SparseColumnLowerFactor & factor,
                                     double modification_parameter = 0.97,
                                     double min_diagonal_ratio = 0.25)
{
  // first copy lower triangle of matrix into factor (Note: assuming A is symmetric of course!)
  factor.resize(matrix.n);
  zero(factor.invdiag); // important: eliminate old values from previous solves!
  factor.value.resize(0);
  factor.rowindex.resize(0);
  zero(factor.adiag);
  for(unsigned int i = 0; (int)i < matrix.n; ++i)
    {
      factor.colstart[i] = (unsigned int)factor.rowindex.size();
      for(unsigned int j = 0; (int)j < matrix.nnz_row(i); ++j)
        {
          if(matrix.col_from_offset(i, j) > (int)i)
            {
              factor.rowindex.push_back(matrix.col_from_offset(i, j));
              factor.value.push_back(matrix.value_from_offset(i, j));
            }
          else if(matrix.col_from_offset(i, j) == (int)i)
            {
              factor.invdiag[i] = factor.adiag[i] = matrix.value_from_offset(i, j);
            }
        }
    }
  factor.colstart[matrix.n] = (unsigned int)factor.rowindex.size();
  // now do the incomplete factorization (figure out numerical values)

  // MATLAB code:
  // L=tril(A);
  // for k=1:size(L,2)
  //   L(k,k)=sqrt(L(k,k));
  //   L(k+1:end,k)=L(k+1:end,k)/L(k,k);
  //   for j=find(L(:,k))'
  //     if j>k
  //       fullupdate=L(:,k)*L(j,k);
  //       incompleteupdate=fullupdate.*(A(:,j)~=0);
  //       missing=sum(fullupdate-incompleteupdate);
  //       L(j:end,j)=L(j:end,j)-incompleteupdate(j:end);
  //       L(j,j)=L(j,j)-omega*missing;
  //     end
  //   end
  // end

  for(unsigned int k = 0; (int)k < matrix.n; ++k)
    {
      if(factor.adiag[k] == 0) continue; // null row/column
      // figure out the final L(k,k) entry
      if(factor.invdiag[k] < min_diagonal_ratio * factor.adiag[k])
        factor.invdiag[k] = 1 / sqrt(factor.adiag[k]); // drop to Gauss-Seidel here if the pivot looks dangerously small
      else
        factor.invdiag[k] = 1 / sqrt(factor.invdiag[k]);
      // finalize the k'th column L(:,k)
      for(unsigned int p = factor.colstart[k]; p < factor.colstart[k + 1]; ++p)
        {
          factor.value[p] *= factor.invdiag[k];
        }
      // incompletely eliminate L(:,k) from future columns, modifying diagonals
      for(unsigned int p = factor.colstart[k]; p < factor.colstart[k + 1]; ++p)
        {
          unsigned int j = factor.rowindex[p]; // work on column j
          double multiplier = factor.value[p];
          double missing = 0;
          unsigned int a = factor.colstart[k];
          // first look for contributions to missing from dropped entries above the diagonal in column j
          unsigned int b = 0;
          while(a < factor.colstart[k + 1] && factor.rowindex[a] < j)
            {
              // look for factor.rowindex[a] in matrix.index[j] starting at b
              while((int)b < matrix.nnz_row(j))
                {
                  if(matrix.col_from_offset(j, b) < (int)factor.rowindex[a])
                    ++b;
                  else if(matrix.col_from_offset(j, b) == (int)factor.rowindex[a])
                    break;
                  else
                    {
                      missing += factor.value[a];
                      break;
                    }
                }
              ++a;
            }
          // adjust the diagonal j,j entry
          if(a < factor.colstart[k + 1] && factor.rowindex[a] == j)
            {
              factor.invdiag[j] -= multiplier * factor.value[a];
            }
          ++a;
          // and now eliminate from the nonzero entries below the diagonal in column j (or add to missing if we can't)
          b = factor.colstart[j];
          while(a < factor.colstart[k + 1] && b < factor.colstart[j + 1])
            {
              if(factor.rowindex[b] < factor.rowindex[a])
                ++b;
              else if(factor.rowindex[b] == factor.rowindex[a])
                {
                  factor.value[b] -= multiplier * factor.value[a];
                  ++a;
                  ++b;
                }
              else
                {
                  missing += factor.value[a];
                  ++a;
                }
            }
          // and if there's anything left to do, add it to missing
          while(a < factor.colstart[k + 1])
            {
              missing += factor.value[a];
              ++a;
            }
          // and do the final diagonal adjustment from the missing entries
          factor.invdiag[j] -= modification_parameter * multiplier * missing;
        }
    }
}

//============================================================================
// Solution routines with lower triangular matrix.

// solve L*result=rhs
inline void
solve_lower(const SparseColumnLowerFactor & factor, const std::vector<double> & rhs, std::vector<double> & result)
{
  assert(factor.n == rhs.size());
  assert(factor.n == result.size());
  result = rhs;
  for(unsigned int i = 0; i < factor.n; ++i)
    {
      result[i] *= factor.invdiag[i];
      for(unsigned int j = factor.colstart[i]; j < factor.colstart[i + 1]; ++j)
        {
          result[factor.rowindex[j]] -= factor.value[j] * result[i];
        }
    }
}

// solve L^T*result=rhs
inline void
solve_lower_transpose_in_place(const SparseColumnLowerFactor & factor, std::vector<double> & x)
{
  assert(factor.n == x.size());
  assert(factor.n > 0);
  unsigned int i = factor.n;
  do
    {
      --i;
      for(unsigned int j = factor.colstart[i]; j < factor.colstart[i + 1]; ++j)
        {
          x[i] -= factor.value[j] * x[factor.rowindex[j]];
        }
      x[i] *= factor.invdiag[i];
    }
  while(i != 0);
}

} // End of namespace bridson

#endif