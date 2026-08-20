#ifndef BRIDSON_SPARSE_MATRIX_H
#define BRIDSON_SPARSE_MATRIX_H

#include "util.h"
#include "vec.h"
#include <algorithm>
#include <climits>
#include <iomanip>
#include <iostream>
#include <vector>

namespace bridson
{


//
// Fixed version of SparseMatrix. This is not a good structure for dynamically
// modifying the matrix, but can be significantly faster for matrix-vector
// multiplies due to better data locality.
//

//
// TODO : Add matmat multiplication.
//

struct FixedSparseMatrix
{
  int n;                      // dimension
  std::vector<double> values; // nonzero values row by row
  std::vector<int> adj;       // corresponding column indices
  std::vector<int> xadj;      // where each row starts in values and adj (and last entry is one past the end, the
                              // number of nonzeros)

  explicit FixedSparseMatrix(int n_ = 0)
  : n(n_)
  , values(0)
  , adj(0)
  , xadj(n_ + 1)
  {
  }

  void clear(void)
  {
    n = 0;
    values.clear();
    adj.clear();
    xadj.clear();
  }

  void resize(int n_)
  {
    n = n_;
    xadj.resize(n + 1);
  }


  void set_adjacency(const std::vector<int> & xadj_in, const std::vector<int> & adj_in)
  {
    xadj = xadj_in;
    adj = adj_in;
    values.resize(adj.size());
    n = (int)xadj.size() - 1;

    for(int i = 0; i < n; ++i)
      {
        assert(xadj[i] >= 0);
        assert(xadj[i] < xadj[i + 1]);
        std::sort(cols_begin(i), cols_end(i));
      }
    assert((int)adj.size() == xadj.back());
    assert(0 == xadj.front());
  }


  //
  // Access specifiers and Editing
  //
  void set_zero() { std::fill(values.begin(), values.end(), double(0)); }

  int nnz_row(const int i) const
  {
    assert(i < n);
    return xadj[i + 1] - xadj[i];
  }

  std::vector<int>::iterator cols_begin(const int i)
  {
    assert(i < n);
    return adj.begin() + xadj[i];
  }

  std::vector<int>::const_iterator cols_cbegin(const int i) const
  {
    assert(i < n);
    return adj.cbegin() + xadj[i];
  }

  std::vector<int>::iterator cols_end(const int i)
  {
    assert(i < n);
    return adj.begin() + xadj[i + 1];
  }

  std::vector<int>::const_iterator cols_cend(const int i) const
  {
    assert(i < n);
    return adj.cbegin() + xadj[i + 1];
  }

  int flatten(const int i, const int j) const
  {
    assert(i < n);
    assert(j < n);
    std::vector<int>::const_iterator flat_iter = std::lower_bound(cols_cbegin(i), cols_cend(i), j);
    assert(flat_iter != adj.end());
    assert(*flat_iter == j);
    return (int)(flat_iter - adj.begin());
  }

  int flatten_from_offset(const int i, const int offset) const
  {
    assert(i < n);
    assert(offset < nnz_row(i));
    return xadj[i] + offset;
  }

  int col_from_offset(const int i, const int offset) const
  {
    assert(i < n);
    assert(offset < nnz_row(i));
    return adj[xadj[i] + offset];
  }

  double value_from_offset(const int i, const int offset) const
  {
    assert(i < n);
    assert(offset < nnz_row(i));
    return values[xadj[i] + offset];
  }

  template<int N, int M>
  void flatten(const Vec<N, int> & row_ids, const Vec<M, int> & col_ids, Vec<N * M, int> & flat_ids) const
  {
    // Can be mad faster, for now ...
    for(int offset_i = 0; offset_i < N; ++offset_i)
      {
        for(int offset_j = 0; offset_j < M; ++offset_j)
          {
            const int fl_gid = flatten(row_ids[offset_i], col_ids[offset_j]);
            const int fl_lid = offset_i * M + offset_j;
            flat_ids[fl_lid] = fl_gid;
          }
      }
  }

  void set_value(const int i, const int j, const double & v) { values[flatten(i, j)] = v; }
  double get_value(const int i, const int j) const { return values[flatten(i, j)]; }

  template<unsigned int N, unsigned int M>
  void set_values(const Vec<N, int> & row_ids, const Vec<M, int> & col_ids, const Vec<N * M, double> & packed_mat)
  {
    for(int offset_i = 0; offset_i < N; ++offset_i)
      {
        for(unsigned int offset_j = 0; offset_j < M; ++offset_j)
          {
            const int fl_gid = flatten(row_ids[offset_i], col_ids[offset_j]);
            const int fl_lid = offset_i * M + offset_j;
            values[fl_gid] = packed_mat[fl_lid];
          }
      }
  }

  template<unsigned int N, unsigned int M>
  void add_values(const Vec<N, int> & row_ids, const Vec<M, int> & col_ids, const Vec<N * M, double> & packed_mat)
  {
    for(unsigned int offset_i = 0; offset_i < N; ++offset_i)
      {
        for(unsigned int offset_j = 0; offset_j < M; ++offset_j)
          {
            const int fl_gid = flatten(row_ids[offset_i], col_ids[offset_j]);
            const int fl_lid = offset_i * M + offset_j;
            values[fl_gid] += packed_mat[fl_lid];
          }
      }
  }

  template<unsigned int N, unsigned int M>
  void get_values(const Vec<N, int> & row_ids, const Vec<M, int> & col_ids, Vec<N * M, double> & packed_mat) const
  {
    for(unsigned int offset_i = 0; offset_i < N; ++offset_i)
      {
        for(unsigned int offset_j = 0; offset_j < M; ++offset_j)
          {
            const int fl_gid = flatten(row_ids[offset_i], col_ids[offset_j]);
            const int fl_lid = offset_i * M + offset_j;
            packed_mat[fl_lid] = values[fl_gid];
          }
      }
  }


  void write_matlab(std::ostream & output, const char * variable_name)
  {
    output << variable_name << "=sparse([";
    for(int i = 0; i < n; ++i)
      {
        for(int j = xadj[i]; j < xadj[i + 1]; ++j)
          {
            output << i + 1 << " ";
          }
      }
    output << "],...\n  [";
    for(int i = 0; i < n; ++i)
      {
        for(int j = xadj[i]; j < xadj[i + 1]; ++j)
          {
            output << adj[j] + 1 << " ";
          }
      }
    output << "],...\n  [";
    for(int i = 0; i < n; ++i)
      {
        for(int j = xadj[i]; j < xadj[i + 1]; ++j)
          {
            output << std::scientific << std::setprecision(16) << values[j] << " ";
          }
      }
    output << "], " << n << ", " << n << ");" << std::endl;
  }
};


// perform result=result+scale*matrix*x
inline void
multiply_scale_and_add(const FixedSparseMatrix & matrix,
                       const std::vector<double> & x,
                       const double scalar,
                       std::vector<double> & result)
{
  assert(matrix.n == (int)x.size());
  result.resize(matrix.n);
  for(unsigned int i = 0; (int)i < matrix.n; ++i)
    {
      for(unsigned int j = matrix.xadj[i]; (int)j < matrix.xadj[i + 1]; ++j)
        {
          result[i] += scalar * matrix.values[j] * x[matrix.adj[j]];
        }
    }
}

// perform result=matrix*x
inline void
multiply(const FixedSparseMatrix & matrix, const std::vector<double> & x, std::vector<double> & result)
{
  assert(matrix.n == (int)x.size());
  result.resize(matrix.n);
  for(int i = 0; i < matrix.n; ++i)
    {
      result[i] = 0;
      for(int flatid = matrix.xadj[i]; flatid < matrix.xadj[i + 1]; ++flatid)
        {
          const int j = matrix.adj[flatid];
          result[i] += matrix.values[flatid] * x[j];
        }
    }
}

// perform result=result+scalar*matrix*x
inline void
multiplyT_scale_and_add(const FixedSparseMatrix & matrix,
                        const std::vector<double> & x,
                        const double scalar,
                        std::vector<double> & result)
{
  assert(matrix.n == (int)x.size());
  result.resize(matrix.n);
  for(int i = 0; i < matrix.n; ++i)
    {
      for(int flatid = matrix.xadj[i]; flatid < matrix.xadj[i + 1]; ++flatid)
        {
          const int j = matrix.adj[flatid];
          result[j] += scalar * matrix.values[flatid] * x[i];
        }
    }
}

// perform result=matrix*x
inline void
multiplyT(const FixedSparseMatrix & matrix, const std::vector<double> & x, std::vector<double> & result)
{
  assert(matrix.n == (int)x.size());
  result.resize(matrix.n);
  std::fill(result.begin(), result.end(), 0);
  multiplyT_scale_and_add(matrix, x, +1, result);
}


} // end of namespace bridson

#endif
