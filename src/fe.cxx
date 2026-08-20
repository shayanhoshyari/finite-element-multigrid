#include "fe.hxx"
#include "cell_type.hxx"

enum
{
  XI = 0,
  ETA
};
enum
{
  X = 0,
  Y
};


FE_segment::FE_segment(const int qorder)
: FE_base()
{
  quadrature().init(CELL_TYPE_EDGE, qorder);
  //
  _JxW.resize(quadrature().n_points());
  _qp.resize(quadrature().n_points());
  //
  _phi.resize(2);
  _phi[0].resize(quadrature().n_points());
  _phi[1].resize(quadrature().n_points());
  //
  _dphi.resize(0);

  // Now start evaluatin
  for(int i = 0; i < (int)quadrature().n_points(); ++i)
    {
      const bridson::Vec3d & xi = quadrature().qp(i);
      _phi[0][i] = (1 - xi(XI)) / 2.;
      _phi[1][i] = (1 + xi(XI)) / 2.;
    }
}

void
FE_segment::reinit(const bridson::Vec2d & v0, const bridson::Vec2d & v1)
{
  // Evaluate the Jacobian (constant in a triangle)
  const double dphi0_dxi(-0.5);
  const double dphi1_dxi(+0.5);

  const bridson::Vec2d dx_dxi = v0 * dphi0_dxi + v1 * dphi1_dxi;

  const double jacobian = bridson::mag(dx_dxi); // not mag2

  // Now start evaluating
  for(int i = 0; i < (int)quadrature().n_points(); ++i)
    {
      _qp[i] = _phi[0][i] * v0 + _phi[1][i] * v1;
      _JxW[i] = jacobian * quadrature().w(i);
   }
}

FE_triangle::FE_triangle(const int qorder)
: FE_base()
{
  quadrature().init(CELL_TYPE_TRI, qorder);
  //
  _JxW.resize(quadrature().n_points());
  _qp.resize(quadrature().n_points());
  //
  _dphi.resize(3);
  _dphi[0].resize(quadrature().n_points());
  _dphi[1].resize(quadrature().n_points());
  _dphi[2].resize(quadrature().n_points());
  //
  _phi.resize(3);
  _phi[0].resize(quadrature().n_points());
  _phi[1].resize(quadrature().n_points());
  _phi[2].resize(quadrature().n_points());

  // Now start evaluating
  for(int i = 0; i < (int)quadrature().n_points(); ++i)
    {
      const bridson::Vec3d & xi = quadrature().qp(i);
      _phi[0][i] = 1 - xi(XI) - xi(ETA);
      _phi[1][i] = xi(XI);
      _phi[2][i] = xi(ETA);
    }
}

void
FE_triangle::reinit(const bridson::Vec2d & v0, const bridson::Vec2d & v1, const bridson::Vec2d & v2)
{
  // Evaluate the Jacobian (constant in a triangle)
  const bridson::Vec2d dphi0_dxi(-1, -1);
  const bridson::Vec2d dphi1_dxi(1, 0);
  const bridson::Vec2d dphi2_dxi(0, 1);

  const bridson::Vec2d dx0_dxi = v0(X) * dphi0_dxi + v1(X) * dphi1_dxi + v2(X) * dphi2_dxi;
  const bridson::Vec2d dx1_dxi = v0(Y) * dphi0_dxi + v1(Y) * dphi1_dxi + v2(Y) * dphi2_dxi;
  const bridson::Vec4d dx_dxi(dx0_dxi(XI), dx0_dxi(ETA), dx1_dxi(XI), dx1_dxi(ETA));

  const double jacobian = dx_dxi(X * 2 + XI) * dx_dxi(Y * 2 + ETA) - dx_dxi(X * 2 + ETA) * dx_dxi(Y * 2 + XI);
  const bridson::Vec4d dxi_dx =
    bridson::Vec4d(dx_dxi(Y * 2 + ETA), -dx_dxi(X * 2 + ETA), -dx_dxi(Y * 2 + XI), dx_dxi(X * 2 + XI)) / jacobian;

  // Now start evaluating
  for(int i = 0; i < (int)quadrature().n_points(); ++i)
    {
      _qp[i] = _phi[0][i] * v0 + _phi[1][i] * v1 + _phi[2][i] * v2;
      _JxW[i] = jacobian * quadrature().w(i);

      _dphi[0][i] = bridson::mat_prodT(dxi_dx, dphi0_dxi);
      _dphi[1][i] = bridson::mat_prodT(dxi_dx, dphi1_dxi);
      _dphi[2][i] = bridson::mat_prodT(dxi_dx, dphi2_dxi);
    }
} // End of reinit
