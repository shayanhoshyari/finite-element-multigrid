#ifndef FE_BASIS_HXX_IS_INCLUDED
#define FE_BASIS_HXX_IS_INCLUDED

#include "bridson/vec.h"
#include "quadrature.hxx"

/**
 * Lagrange shape functions
 */
class FE_base
{
public:
  const std::vector<bridson::Vec2d> & qp() const { return _qp; }
  const std::vector<double> & JxW() const { return _JxW; }
  // Index 0 dof, Index 1 qp
  const std::vector<std::vector<double>> & phi() const { return _phi; }
  const std::vector<std::vector<bridson::Vec2d>> & dphi() const { return _dphi; }

protected:
  FE_base() = default;
  std::vector<bridson::Vec2d> _qp;
  std::vector<double> _JxW;
  std::vector<std::vector<double> > _phi;
  std::vector<std::vector<bridson::Vec2d> > _dphi;
  QGenerator _quadrature;

  //std::vector<bridson::Vec2d> & qp() { return _qp; }
  //std::vector<double> & JxW() { return _JxW; }
  //std::vector<std::vector<double> > & phi() { return _phi; }
  //std::vector<std::vector<bridson::Vec2d> > & dphi() { return _dphi; }
  QGenerator & quadrature() { return _quadrature; }
};

class FE_segment : public FE_base
{
public:
  FE_segment(const int qorder);
  void reinit(const bridson::Vec2d & v0, const bridson::Vec2d & v1);
private:
};

class FE_triangle : public FE_base
{
public:
  FE_triangle(const int qorder);
  void reinit(const bridson::Vec2d & v0, const bridson::Vec2d & v1, const bridson::Vec2d & v2);

};


#endif