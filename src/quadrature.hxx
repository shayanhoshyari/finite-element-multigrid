#ifndef QUADRATURE_HXX_IS_INCLUDED
#define QUADRATURE_HXX_IS_INCLUDED

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "cell_type.hxx"
#include "bridson/vec.h"

/* -----------------------------------------------------------------------
 *  Point
 * ----------------------------------------------------------------------- */


class QGenerator
{
  typedef bridson::Vec3d Point;
  std::vector<Point> _points;
  std::vector<double> _weights;

  CellType _type;
  unsigned int _order;

  // scale in the 1D case
  // assumes old range is [-1 1]
  void scale(const double new_range[]);

  // quadrature rules for triangle
  void dunavant_rule2(const double * wts,
                      const double * a,
                      const double * b,
                      const unsigned int * permutation_ids,
                      const unsigned int n_wts);

  // quadrature rule for quad
  void tensor_product_quad();

  // quadrature rule for tet
  void conical_rule_tet();

  // quadrature rule for prism
  void tensor_product_prism();

  // quadrature rule for pyramid
  void conical_rule_pyramid();

  // quadrature rule for hex
  void tensor_product_hex();



public:
  unsigned int n_points();
  const Point & qp(const int i) const;
  double w(const int i) const;


  void init(CellType type, const unsigned int order);
  void init_1d_jacobi(const int a);
  void init_1d_gauss();
  void init_2d_gauss();
  void init_3d_gauss();
};


#endif /* QUADRATURE_HXX_IS_INCLUDED */
