#ifndef NUMERIC_ARRAY_IS_INCLUDED
#define NUMERIC_ARRAY_IS_INCLUDED

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include "bridson/blas_wrappers.h"
#include "bridson/vec.h"


class Numeric_array : public std::vector<double>
{
private:
  bool _is_init;
  unsigned _fixed_size;
  typedef std::vector<double> Parent_type;
  void verify() const { assert(_fixed_size == Parent_type::size() && "Should not change size"); }

public:
  Numeric_array(const unsigned size)
  : Parent_type(size)
  , _is_init(true)
  , _fixed_size(size)
  {
  }

  Numeric_array()
  : Parent_type()
  , _is_init(false)
  , _fixed_size(0)
  {
  }

  // This must be explicit to prevent human and compiler confusion
  explicit Numeric_array(const Parent_type & p)
  : Parent_type(p)
  , _is_init(true)
  , _fixed_size((unsigned)p.size())
  {
  }

  void init(const unsigned size)
  {
    assert(!_is_init);
    Parent_type::resize(size);
    _fixed_size = size;
    _is_init = true;
  }

  void init(const Parent_type & p)
  {
    assert(!_is_init);
    Parent_type::resize(p.size());
    _fixed_size = (unsigned)p.size();
    std::copy(p.begin(), p.end(), Parent_type::begin());
    _is_init = true;
  }

  void init_swap(Parent_type & p)
  {
    assert(!_is_init);
    _fixed_size = (unsigned)p.size();
    Parent_type::swap(p);
    _is_init = true;
  }

  void swap(Parent_type & p)
  {
    assert(_fixed_size == (unsigned)p.size());
    Parent_type::swap(p);
  }

  void add(const Parent_type & o)
  {
    verify();
    add_scaled(+1, o);
  }

  void subtract(const Parent_type & o)
  {
    verify();
    add_scaled(-1, o);
  }

  void add_scaled(const double a, const Parent_type & o)
  {
    verify();
    bridson::BLAS::add_scaled(a, o, *this);
  }

  void set_equal(const double val)
  {
    verify();
    std::fill(this->begin(), this->end(), val);
  }

  void set_zero()
  {
    verify();
    this->set_equal(0);
  }

  void set_equal(const Parent_type & o)
  {
    verify();
    assert(this->size() == o.size());
    std::copy(o.begin(), o.end(), this->begin());
  }

  int arg_max() const
  {
    verify();
    return bridson::BLAS::index_abs_max(*this);
  }

  double dot(const Parent_type & o) const
  {
    verify();
    return bridson::BLAS::dot(*this, *this);
  }

  double norm_l2_squared() const
  {
    verify();
    return dot(*this);
  }

  double norm_l2() const
  {
    verify();
    return sqrt(norm_l2_squared());
  }

  double norm_linf() const
  {
    verify();
    return bridson::BLAS::abs_max(*this);
  }

  double norm_l1() const
  {
    verify();
    return bridson::BLAS::sum_abs(*this);
  }

  template<unsigned int N>
  void add_values(const bridson::Vec<N, int> & ids, const bridson::Vec<N, double> & vals)
  {
    verify();
    for(unsigned int i = 0; i < N; ++i)
      {
        assert(ids[i] < (int)size());
        this->at(ids[i]) += vals[i];
      }
  }

  std::string to_string() const
  {
    verify();
    std::stringstream ss;
    Parent_type::const_iterator it = this->cbegin();
    for(; it != this->cend(); ++it)
      {
        ss << *it << "\n";
      }
    return ss.str();
  }

  void write_matlab(std::ostream & output, const char * variable_name)
  {
    verify();
    output << variable_name << " = [";
    for(Parent_type::iterator it = this->begin(); it != this->end(); ++it)
      {
        output << std::scientific << std::setprecision(16) << *it << ", ";
      }
    output << "]';" << std::endl;
  }
};

#endif /*NUMERIC_ARRAY_IS_INCLUDED*/
