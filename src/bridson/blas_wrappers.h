#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

#include <cmath>
#include <vector>

#ifdef USE_CBLAS
#include <cblas.h>
#endif

namespace bridson
{
namespace BLAS
{

// dot products ==============================================================

inline double
dot(const std::vector<double> & x, const std::vector<double> & y)
{
#ifdef USE_CBLAS
   return cblas_ddot((int)x.size(), x.data(), 1, y.data(), 1);
#else
  double sum = 0;
  for(unsigned int i = 0; i < x.size(); ++i) sum += x[i] * y[i];
  return sum;
#endif
}

// inf-norm (maximum absolute value: index of max returned) ==================

inline int
index_abs_max(const std::vector<double> & x)
{
#ifdef USE_CBLAS
  return cblas_idamax((int)x.size(), x.data(), 1);
#else
  int maxind = 0;
  double maxvalue = 0;
  for(unsigned int i = 0; i < x.size(); ++i)
    {
      if(std::abs((double)x[i]) > maxvalue)
        {
          maxvalue = std::abs((double)x[i]);
          maxind = i;
        }
    }
  return maxind;
#endif
}

// inf-norm (maximum absolute value) =========================================
// technically not part of BLAS, but useful
inline double
abs_max(const std::vector<double> & x)
{
  return std::abs(x[index_abs_max(x)]);
}

// saxpy (y=alpha*x+y) =======================================================
inline void
add_scaled(double alpha, const std::vector<double> & x, std::vector<double> & y)
{
#ifdef USE_CBLAS
  cblas_daxpy((int)x.size(), alpha, x.data(), 1, y.data(), 1);
#else
  for(unsigned int i = 0; i < x.size(); ++i)
  {
    y[i] = alpha*x[i]+y[i];
  }
#endif
}

inline double
sum_abs(const std::vector<double> & x)
{
#ifdef USE_CBLAS
  return cblas_dasum((int)x.size(), x.data(), 1);
#else
  double sum = 0;
  for(unsigned int i = 0; i < x.size(); ++i) sum += std::abs(x[i]);
  return sum;
#endif
}
}
}
#endif
