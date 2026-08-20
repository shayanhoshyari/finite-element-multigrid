/*
 * The following code is taken from the libMesh library, with slight
 * modifications. libMesh/src/quadrature/
 */

#include <cassert>

#include "quadrature.hxx"


void
QGenerator::init(CellType type, const unsigned int order)
{
  _type = type;
  _order = order;
  if(cell_dim(type) == 1) init_1d_gauss();
  else if(cell_dim(type) == 2) init_2d_gauss();
  else if(cell_dim(type) == 3) init_3d_gauss();
  else assert(0 && "Dimension not supported");
}

unsigned int
QGenerator::n_points()
{
  return _points.size();
}

const QGenerator::Point &
QGenerator::qp(const int i) const
{
  return _points[i];
}

double
QGenerator::w(const int i) const
{
  return _weights[i];
}

/* -----------------------------------------------------------------------
 *   1D rules
 * ----------------------------------------------------------------------- */
void
QGenerator::scale(const double new_range[])
{
  // Make sure we are in 1D
  assert(cell_dim(_type) == 1);

  double h_new = new_range[1] - new_range[0], h_old = 2;

  // Make sure that we have sane ranges
  assert(h_new > 0.);

  // Make sure there are some points
  assert(_points.size() > 0);

  // Compute the scale factor
  double scfact = h_new / h_old;

  // We're mapping from old_range -> new_range
  for(unsigned int i = 0; i < _points.size(); i++)
    {
      _points[i](0) = new_range[0] + (_points[i](0) - (-1)) * scfact;

      // Scale the weights
      _weights[i] *= scfact;
    }
}

void
QGenerator::init_1d_gauss()
{

  assert(_type == CELL_TYPE_EDGE);

  switch(_order)
    {
      case 0:
      case 1:
        {
          _points.resize(1);
          _weights.resize(1);

          _points[0](0) = 0.;

          _weights[0] = 2.;

          return;
        }
      case 2:
      case 3:
        {
          _points.resize(2);
          _weights.resize(2);

          _points[0](0) = -5.7735026918962576450914878050196e-01L; // -sqrt(3)/3
          _points[1] = -_points[0];

          _weights[0] = 1.;
          _weights[1] = _weights[0];

          return;
        }
      case 4:
      case 5:
        {
          _points.resize(3);
          _weights.resize(3);

          _points[0](0) = -7.7459666924148337703585307995648e-01L;
          _points[1](0) = 0.;
          _points[2] = -_points[0];

          _weights[0] = 5.5555555555555555555555555555556e-01L;
          _weights[1] = 8.8888888888888888888888888888889e-01L;
          _weights[2] = _weights[0];

          return;
        }
      case 6:
      case 7:
        {
          _points.resize(4);
          _weights.resize(4);

          _points[0](0) = -8.6113631159405257522394648889281e-01L;
          _points[1](0) = -3.3998104358485626480266575910324e-01L;
          _points[2] = -_points[1];
          _points[3] = -_points[0];

          _weights[0] = 3.4785484513745385737306394922200e-01L;
          _weights[1] = 6.5214515486254614262693605077800e-01L;
          _weights[2] = _weights[1];
          _weights[3] = _weights[0];

          return;
        }
      case 8:
      case 9:
        {
          _points.resize(5);
          _weights.resize(5);

          _points[0](0) = -9.0617984593866399279762687829939e-01L;
          _points[1](0) = -5.3846931010568309103631442070021e-01L;
          _points[2](0) = 0.;
          _points[3] = -_points[1];
          _points[4] = -_points[0];

          _weights[0] = 2.3692688505618908751426404071992e-01L;
          _weights[1] = 4.7862867049936646804129151483564e-01L;
          _weights[2] = 5.6888888888888888888888888888889e-01L;
          _weights[3] = _weights[1];
          _weights[4] = _weights[0];

          return;
        }
      case 10:
      case 11:
        {
          _points.resize(6);
          _weights.resize(6);

          _points[0](0) = -9.3246951420315202781230155449399e-01L;
          _points[1](0) = -6.6120938646626451366139959501991e-01L;
          _points[2](0) = -2.3861918608319690863050172168071e-01L;
          _points[3] = -_points[2];
          _points[4] = -_points[1];
          _points[5] = -_points[0];

          _weights[0] = 1.7132449237917034504029614217273e-01L;
          _weights[1] = 3.6076157304813860756983351383772e-01L;
          _weights[2] = 4.6791393457269104738987034398955e-01L;
          _weights[3] = _weights[2];
          _weights[4] = _weights[1];
          _weights[5] = _weights[0];

          return;
        }

      default: assert("Order not supported for Gauss" && 0);
    }
}

void
QGenerator::init_1d_jacobi(const int a)
{
  switch(a)
    {
      case 1:
        {
          switch(_order)
            {
              case 0:
              case 1:
                {
                  _points.resize(1);
                  _weights.resize(1);

                  _points[0](0) = 1. / 3.;

                  _weights[0] = 0.5;

                  return;
                }
              case 2:
              case 3:
                {
                  _points.resize(2);
                  _weights.resize(2);

                  _points[0](0) = 1.5505102572168219018027159252941e-01L;
                  _points[1](0) = 6.4494897427831780981972840747059e-01L;

                  _weights[0] = 3.1804138174397716939436900207516e-01L;
                  _weights[1] = 1.8195861825602283060563099792484e-01L;

                  return;
                }
              case 4:
              case 5:
                {
                  _points.resize(3);
                  _weights.resize(3);

                  _points[0](0) = 8.8587959512703947395546143769456e-02L;
                  _points[1](0) = 4.0946686444073471086492625206883e-01L;
                  _points[2](0) = 7.8765946176084705602524188987600e-01L;

                  _weights[0] = 2.0093191373895963077219813326460e-01L;
                  _weights[1] = 2.2924110635958624669392059455632e-01L;
                  _weights[2] = 6.9826979901454122533881272179077e-02L;

                  return;
                }
              case 6:
              case 7:
                {
                  _points.resize(4);
                  _weights.resize(4);

                  _points[0](0) = 5.7104196114517682193121192554116e-02L;
                  _points[1](0) = 2.7684301363812382768004599768563e-01L;
                  _points[2](0) = 5.8359043236891682005669766866292e-01L;
                  _points[3](0) = 8.6024013565621944784791291887512e-01L;

                  _weights[0] = 1.3550691343148811620826417407794e-01L;
                  _weights[1] = 2.0346456801027136079140447593585e-01L;
                  _weights[2] = 1.2984754760823244082645620288963e-01L;
                  _weights[3] = 3.1180970950008082173875147096575e-02L;

                  return;
                }
              case 8:
              case 9:
                {
                  _points.resize(5);
                  _weights.resize(5);

                  _points[0](0) = 3.9809857051468742340806690093333e-02L;
                  _points[1](0) = 1.9801341787360817253579213679530e-01L;
                  _points[2](0) = 4.3797481024738614400501252000523e-01L;
                  _points[3](0) = 6.9546427335363609451461482372117e-01L;
                  _points[4](0) = 9.0146491420117357387650110211225e-01L;

                  _weights[0] = 9.6781590226651679274360971636151e-02L;
                  _weights[1] = 1.6717463809436956549167562309770e-01L;
                  _weights[2] = 1.4638698708466980869803786935596e-01L;
                  _weights[3] = 7.3908870072616670350633219341704e-02L;
                  _weights[4] = 1.5747914521692276185292316568488e-02L;

                  return;
                }
              case 10:
              case 11:
                {

                  _points.resize(6);
                  _weights.resize(6);

                  _points[0](0) = 2.9316427159784891972050276913165e-02L;
                  _points[1](0) = 1.4807859966848429184997685249598e-01L;
                  _points[2](0) = 3.3698469028115429909705297208078e-01L;
                  _points[3](0) = 5.5867151877155013208139334180552e-01L;
                  _points[4](0) = 7.6923386203005450091688336011565e-01L;
                  _points[5](0) = 9.2694567131974111485187396581968e-01L;

                  _weights[0] = 7.2310330725508683655454326124832e-02L;
                  _weights[1] = 1.3554249723151861684069039663805e-01L;
                  _weights[2] = 1.4079255378819892811907683907093e-01L;
                  _weights[3] = 9.8661150890655264120584510548360e-02L;
                  _weights[4] = 4.3955165550508975508176624305430e-02L;
                  _weights[5] = 8.7383018136095317560173033123989e-03L;

                  return;
                }

              default: assert("Order not supported for jacobi a=1" && 0);
            } // end switch(_order)

        } // end case(a=1)


      case 2:
        {
          switch(_order)
            {
              case 0:
              case 1:
                {
                  _points.resize(1);
                  _weights.resize(1);

                  _points[0](0) = 0.25;

                  _weights[0] = 1. / 3.;

                  return;
                }
              case 2:
              case 3:
                {
                  _points.resize(2);
                  _weights.resize(2);

                  _points[0](0) = 1.2251482265544137786674043037115e-01L;
                  _points[1](0) = 5.4415184401122528879992623629551e-01L;

                  _weights[0] = 2.3254745125350790274997694884235e-01L;
                  _weights[1] = 1.0078588207982543058335638449099e-01L;

                  return;
                }
              case 4:
              case 5:
                {
                  _points.resize(3);
                  _weights.resize(3);

                  _points[0](0) = 7.2994024073149732155837979012003e-02L;
                  _points[1](0) = 3.4700376603835188472176354340395e-01L;
                  _points[2](0) = 7.0500220988849838312239847758405e-01L;

                  _weights[0] = 1.5713636106488661331834482221327e-01L;
                  _weights[1] = 1.4624626925986602200351202036424e-01L;
                  _weights[2] = 2.9950703008580698011476490755827e-02L;

                  return;
                }
              case 6:
              case 7:
                {
                  _points.resize(4);
                  _weights.resize(4);

                  _points[0](0) = 4.8500549446997329297067257098986e-02L;
                  _points[1](0) = 2.3860073755186230505898141272470e-01L;
                  _points[2](0) = 5.1704729510436750234057336938307e-01L;
                  _points[3](0) = 7.9585141789677286330337796079324e-01L;

                  _weights[0] = 1.1088841561127798368323131746895e-01L;
                  _weights[1] = 1.4345878979921420904832801427594e-01L;
                  _weights[2] = 6.8633887172923075317376345041811e-02L;
                  _weights[3] = 1.0352240749918065284397656546639e-02L;

                  return;
                }
              case 8:
              case 9:
                {
                  _points.resize(5);
                  _weights.resize(5);

                  _points[0](0) = 3.4578939918215091524457428631527e-02L;
                  _points[1](0) = 1.7348032077169572310459241798618e-01L;
                  _points[2](0) = 3.8988638706551932824089541038499e-01L;
                  _points[3](0) = 6.3433347263088677234716388892062e-01L;
                  _points[4](0) = 8.5105421294701641811622418741001e-01L;

                  _weights[0] = 8.1764784285770917904880732922352e-02L;
                  _weights[1] = 1.2619896189991148802883293516467e-01L;
                  _weights[2] = 8.9200161221590000186254493070384e-02L;
                  _weights[3] = 3.2055600722961919254748930556633e-02L;
                  _weights[4] = 4.1138252030990079586162416192983e-03L;

                  return;
                }
              case 10:
              case 11:
                {
                  _points.resize(6);
                  _weights.resize(6);

                  _points[0](0) = 2.5904555093667192754643606997235e-02L;
                  _points[1](0) = 1.3156394165798513398691085074097e-01L;
                  _points[2](0) = 3.0243691802289123274990557791855e-01L;
                  _points[3](0) = 5.0903641316475208401103990516772e-01L;
                  _points[4](0) = 7.1568112731171391876766262459361e-01L;
                  _points[5](0) = 8.8680561617756186630126600601049e-01L;

                  _weights[0] = 6.2538702726580937878526556468332e-02L;
                  _weights[1] = 1.0737649973678063260575568234795e-01L;
                  _weights[2] = 9.4577186748541203568292051720052e-02L;
                  _weights[3] = 5.1289571129616210220129325076919e-02L;
                  _weights[4] = 1.5720297184945051327851262130020e-02L;
                  _weights[5] = 1.8310758068692977327784555900562e-03L;

                  return;
                }

              default: assert("Order not supported for jacobi a=2" && 0);

            } // end switch(_order)
        }     // end case(a=2)
    }         // end switch(a)
}

/* -----------------------------------------------------------------------
 *   2D rules
 * ----------------------------------------------------------------------- */

// A number of different rules for triangles can be described by
// permutations of the following types of points:
// I:   "1"-permutation, (1/3,1/3)  (single point only)
// II:   3-permutation, (a,a,1-2a)
// III:  6-permutation, (a,b,1-a-b)
// The weights for a given set of permutations are all the same.
void
QGenerator::dunavant_rule2(const double * wts,
                           const double * a,
                           const double * b,
                           const unsigned int * permutation_ids,
                           unsigned int n_wts)
{
  // Figure out how many total points by summing up the entries
  // in the permutation_ids array, and resize the _points and _weights
  // vectors appropriately.
  unsigned int total_pts = 0;
  for(unsigned int p = 0; p < n_wts; ++p) total_pts += permutation_ids[p];

  // Resize point and weight vectors appropriately.
  _points.resize(total_pts);
  _weights.resize(total_pts);

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset = 0;

  for(unsigned int p = 0; p < n_wts; ++p)
    {
      switch(permutation_ids[p])
        {
          case 1:
            {
              // The point has only a single permutation (the centroid!)
              // So we don't even need to look in the a or b arrays.
              _points[offset + 0] = Point(1.0L / 3.0L, 1.0L / 3.0L);
              _weights[offset + 0] = wts[p];

              offset += 1;
              break;
            }


          case 3:
            {
              // For this type of rule, don't need to look in the b array.
              _points[offset + 0] = Point(a[p], a[p]);             // (a,a)
              _points[offset + 1] = Point(a[p], 1.L - 2.L * a[p]); // (a,1-2a)
              _points[offset + 2] = Point(1.L - 2.L * a[p], a[p]); // (1-2a,a)

              for(unsigned int j = 0; j < 3; ++j) _weights[offset + j] = wts[p];

              offset += 3;
              break;
            }

          case 6:
            {
              // This type of point uses all 3 arrays...
              _points[offset + 0] = Point(a[p], b[p]);
              _points[offset + 1] = Point(b[p], a[p]);
              _points[offset + 2] = Point(a[p], 1.L - a[p] - b[p]);
              _points[offset + 3] = Point(1.L - a[p] - b[p], a[p]);
              _points[offset + 4] = Point(b[p], 1.L - a[p] - b[p]);
              _points[offset + 5] = Point(1.L - a[p] - b[p], b[p]);

              for(unsigned int j = 0; j < 6; ++j) _weights[offset + j] = wts[p];

              offset += 6;
              break;
            }

          default: assert("Unknown permutation id" && 0);
        }
    }
}

void
QGenerator::tensor_product_quad()
{
  QGenerator q1D;
  q1D.init(CELL_TYPE_EDGE, _order);
  q1D.init_1d_gauss();

  const unsigned int np = q1D.n_points();

  _points.resize(np * np);

  _weights.resize(np * np);

  unsigned int q = 0;

  for(unsigned int j = 0; j < np; j++)
    for(unsigned int i = 0; i < np; i++)
      {
        _points[q](0) = q1D.qp(i)(0);
        _points[q](1) = q1D.qp(j)(0);

        _weights[q] = q1D.w(i) * q1D.w(j);

        q++;
      }
}

void
QGenerator::init_2d_gauss()
{
  switch(_type)
    {

        //---------------------------------------------
        // Quadrilateral quadrature rules
      case CELL_TYPE_QUAD:
        {
          tensor_product_quad();
          return;
        }


        //---------------------------------------------
        // Triangle quadrature rules
      case CELL_TYPE_TRI:
        {
          switch(_order)
            {
              case 0:
              case 1:
                {
                  // Exact for linears
                  _points.resize(1);
                  _weights.resize(1);

                  _points[0](0) = double(1) / 3;
                  _points[0](1) = double(1) / 3;

                  _weights[0] = 0.5;

                  return;
                }
              case 2:
                {
                  // Exact for quadratics
                  _points.resize(3);
                  _weights.resize(3);

                  _points[0](0) = double(2) / 3;
                  _points[0](1) = double(1) / 6;

                  _points[1](0) = double(1) / 6;
                  _points[1](1) = double(2) / 3;

                  _points[2](0) = double(1) / 6;
                  _points[2](1) = double(1) / 6;


                  _weights[0] = double(1) / 6;
                  _weights[1] = double(1) / 6;
                  _weights[2] = double(1) / 6;

                  return;
                }
              case 3:
                {
                  // Exact for cubics
                  _points.resize(4);
                  _weights.resize(4);

                  _points[0](0) = 1.5505102572168219018027159252941e-01L;
                  _points[0](1) = 1.7855872826361642311703513337422e-01L;
                  _points[1](0) = 6.4494897427831780981972840747059e-01L;
                  _points[1](1) = 7.5031110222608118177475598324603e-02L;
                  _points[2](0) = 1.5505102572168219018027159252941e-01L;
                  _points[2](1) = 6.6639024601470138670269327409637e-01L;
                  _points[3](0) = 6.4494897427831780981972840747059e-01L;
                  _points[3](1) = 2.8001991549907407200279599420481e-01L;

                  _weights[0] = 1.5902069087198858469718450103758e-01L;
                  _weights[1] = 9.0979309128011415302815498962418e-02L;
                  _weights[2] = 1.5902069087198858469718450103758e-01L;
                  _weights[3] = 9.0979309128011415302815498962418e-02L;

                  return;
                }

              case 4:
                {
                  const unsigned int n_wts = 2;
                  const double wts[n_wts] = {1.1169079483900573284750350421656140e-01L,
                                             5.4975871827660933819163162450105264e-02L};

                  const double a[n_wts] = {4.4594849091596488631832925388305199e-01L,
                                           9.1576213509770743459571463402201508e-02L};

                  const double b[n_wts] = {0., 0.}; // not used
                  const unsigned int permutation_ids[n_wts] = {3, 3};

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts); // 6 total points

                  return;
                }

              case 5:
                {
                  const unsigned int n_wts = 3;
                  const double wts[n_wts] = {double(9) / 80,
                                             static_cast<double>(double(31) / 480 + std::sqrt(15.0L) / 2400),
                                             static_cast<double>(double(31) / 480 - std::sqrt(15.0L) / 2400)};

                  const double a[n_wts] = {0., // 'a' parameter not used for origin permutation
                                           static_cast<double>(double(2) / 7 + std::sqrt(15.0L) / 21),
                                           static_cast<double>(double(2) / 7 - std::sqrt(15.0L) / 21)};

                  const double b[n_wts] = {0., 0., 0.}; // not used
                  const unsigned int permutation_ids[n_wts] = {1, 3, 3};

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts); // 7 total points

                  return;
                }



              // A degree 6 rule with 12 points.  This rule can be found in many places
              // including:
              //
              // J.N. Lyness and D. Jespersen, Moderate degree symmetric
              // quadrature rules for the triangle, J. Inst. Math. Appl.  15 (1975),
              // 19--32.
              //
              // We used the code in:
              // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
              // on triangles and tetrahedra"  Journal of Computational Mathematics,
              // v. 27, no. 1, 2009, pp. 89-96.
              // to generate additional precision.
              //
              // Note that the following 7th-order Ro3-invariant rule also has only 12 points,
              // which technically makes it the superior rule.  This one is here for completeness.
              case 6:
                {
                  const unsigned int n_wts = 3;
                  const double wts[n_wts] = {5.8393137863189683012644805692789721e-02L,
                                             2.5422453185103408460468404553434492e-02L,
                                             4.1425537809186787596776728210221227e-02L};

                  const double a[n_wts] = {2.4928674517091042129163855310701908e-01L,
                                           6.3089014491502228340331602870819157e-02L,
                                           3.1035245103378440541660773395655215e-01L};

                  const double b[n_wts] = {0., 0., 6.3650249912139864723014259441204970e-01L};

                  const unsigned int permutation_ids[n_wts] = {3, 3, 6}; // 12 total points

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts);

                  return;
                }


              // A degree 7 rule with 12 points.  This rule can be found in:
              //
              // K. Gatermann, The construction of symmetric cubature
              // formulas for the square and the triangle, Computing 40
              // (1988), 229--240.
              //
              // This rule, which is provably minimal in the number of
              // integration points, is said to be 'Ro3 invariant' which
              // means that a given set of barycentric coordinates
              // (z1,z2,z3) implies the quadrature points (z1,z2),
              // (z3,z1), (z2,z3) which are formed by taking the first
              // two entries in cyclic permutations of the barycentric
              // point.  Barycentric coordinates are related in the
              // sense that: z3 = 1 - z1 - z2.
              //
              // The 12-point sixth-order rule for triangles given in
              // Flaherty's (http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps)
              // lecture notes has been removed in favor of this rule
              // which is higher-order (for the same number of
              // quadrature points) and has a few more digits of
              // precision in the points and weights.  Some 10-point
              // degree 6 rules exist for the triangle but they have
              // quadrature points outside the region of integration.
              case 7:
                {
                  _points.resize(12);
                  _weights.resize(12);

                  const unsigned int nrows = 4;

                  // In each of the rows below, the first two entries are (z1, z2) which imply
                  // z3.  The third entry is the weight for each of the points in the cyclic permutation.
                  const double rule_data[nrows][3] = {
                    {6.2382265094402118e-02, 6.7517867073916085e-02, 2.6517028157436251e-02}, // group A
                    {5.5225456656926611e-02, 3.2150249385198182e-01, 4.3881408714446055e-02}, // group B
                    {3.4324302945097146e-02, 6.6094919618673565e-01, 2.8775042784981585e-02}, // group C
                    {5.1584233435359177e-01, 2.7771616697639178e-01, 6.7493187009802774e-02}  // group D
                  };

                  for(unsigned int i = 0, offset = 0; i < nrows; ++i)
                    {
                      _points[offset + 0] = Point(rule_data[i][0], rule_data[i][1]);                        // (z1,z2)
                      _points[offset + 1] = Point(1. - rule_data[i][0] - rule_data[i][1], rule_data[i][0]); // (z3,z1)
                      _points[offset + 2] = Point(rule_data[i][1], 1. - rule_data[i][0] - rule_data[i][1]); // (z2,z3)

                      // All these points get the same weight
                      _weights[offset + 0] = rule_data[i][2];
                      _weights[offset + 1] = rule_data[i][2];
                      _weights[offset + 2] = rule_data[i][2];

                      // Increment offset
                      offset += 3;
                    }

                  return;
                }



              // Another Dunavant rule.  This one has all positive weights.  This rule has
              // 16 points while a comparable conical product rule would have 5*5=25.
              //
              // It was copied 23rd June 2008 from:
              // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
              //
              // Additional precision obtained from the code in:
              // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
              // on triangles and tetrahedra"  Journal of Computational Mathematics,
              // v. 27, no. 1, 2009, pp. 89-96.
              case 8:
                {
                  const unsigned int n_wts = 5;
                  const double wts[n_wts] = {7.2157803838893584125545555244532310e-02L,
                                             4.7545817133642312396948052194292159e-02L,
                                             5.1608685267359125140895775146064515e-02L,
                                             1.6229248811599040155462964170890299e-02L,
                                             1.3615157087217497132422345036954462e-02L};

                  const double a[n_wts] = {
                    0.0, // 'a' parameter not used for origin permutation
                    4.5929258829272315602881551449416932e-01L,
                    1.7056930775176020662229350149146450e-01L,
                    5.0547228317030975458423550596598947e-02L,
                    2.6311282963463811342178578628464359e-01L,
                  };

                  const double b[n_wts] = {0., 0., 0., 0., 7.2849239295540428124100037917606196e-01L};

                  const unsigned int permutation_ids[n_wts] = {1, 3, 3, 3, 6}; // 16 total points

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts);

                  return;
                }



              // Another Dunavant rule.  This one has all positive weights.  This rule has 19
              // points. The comparable conical product rule would have 25.
              // It was copied 23rd June 2008 from:
              // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
              //
              // Additional precision obtained from the code in:
              // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
              // on triangles and tetrahedra"  Journal of Computational Mathematics,
              // v. 27, no. 1, 2009, pp. 89-96.
              case 9:
                {
                  const unsigned int n_wts = 6;
                  const double wts[n_wts] = {4.8567898141399416909620991253644315e-02L,
                                             1.5667350113569535268427415643604658e-02L,
                                             1.2788837829349015630839399279499912e-02L,
                                             3.8913770502387139658369678149701978e-02L,
                                             3.9823869463605126516445887132022637e-02L,
                                             2.1641769688644688644688644688644689e-02L};

                  const double a[n_wts] = {0.0, // 'a' parameter not used for origin permutation
                                           4.8968251919873762778370692483619280e-01L,
                                           4.4729513394452709865106589966276365e-02L,
                                           4.3708959149293663726993036443535497e-01L,
                                           1.8820353561903273024096128046733557e-01L,
                                           2.2196298916076569567510252769319107e-01L};

                  const double b[n_wts] = {0., 0., 0., 0., 0., 7.4119859878449802069007987352342383e-01L};

                  const unsigned int permutation_ids[n_wts] = {1, 3, 3, 3, 3, 6}; // 19 total points

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts);

                  return;
                }


              // Another Dunavant rule with all positive weights.  This rule has 25
              // points. The comparable conical product rule would have 36.
              // It was copied 23rd June 2008 from:
              // http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.f90
              //
              // Additional precision obtained from the code in:
              // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
              // on triangles and tetrahedra"  Journal of Computational Mathematics,
              // v. 27, no. 1, 2009, pp. 89-96.
              case 10:
                {
                  const unsigned int n_wts = 6;
                  const double wts[n_wts] = {4.5408995191376790047643297550014267e-02L,
                                             1.8362978878233352358503035945683300e-02L,
                                             2.2660529717763967391302822369298659e-02L,
                                             3.6378958422710054302157588309680344e-02L,
                                             1.4163621265528742418368530791049552e-02L,
                                             4.7108334818664117299637354834434138e-03L};

                  const double a[n_wts] = {0.0, // 'a' parameter not used for origin permutation
                                           4.8557763338365737736750753220812615e-01L,
                                           1.0948157548503705479545863134052284e-01L,
                                           3.0793983876412095016515502293063162e-01L,
                                           2.4667256063990269391727646541117681e-01L,
                                           6.6803251012200265773540212762024737e-02L};

                  const double b[n_wts] = {0.,
                                           0.,
                                           0.,
                                           5.5035294182099909507816172659300821e-01L,
                                           7.2832390459741092000873505358107866e-01L,
                                           9.2365593358750027664630697761508843e-01L};

                  const unsigned int permutation_ids[n_wts] = {1, 3, 3, 6, 6, 6}; // 25 total points

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts);

                  return;
                }


              // Dunavant's 11th-order rule contains points outside the region of
              // integration, and is thus unacceptable for our FEM calculations.
              //
              // This 30-point, 11th-order rule was obtained by me [JWP] using the code in
              //
              // Additional precision obtained from the code in:
              // L. Zhang, T. Cui, and H. Liu. "A set of symmetric quadrature rules
              // on triangles and tetrahedra"  Journal of Computational Mathematics,
              // v. 27, no. 1, 2009, pp. 89-96.
              //
              // Note: the 28-point 11th-order rule obtained by Zhang in the paper above
              // does not appear to be unique.  It is a solution in the sense that it
              // minimizes the error in the least-squares minimization problem, but
              // it involves too many unknowns and the Jacobian is therefore singular
              // when attempting to improve the solution via Newton's method.
              case 11:
                {
                  const unsigned int n_wts = 6;
                  const double wts[n_wts] = {3.6089021198604635216985338480426484e-02L,
                                             2.1607717807680420303346736867931050e-02L,
                                             3.1144524293927978774861144478241807e-03L,
                                             2.9086855161081509446654185084988077e-02L,
                                             8.4879241614917017182977532679947624e-03L,
                                             1.3795732078224796530729242858347546e-02L};

                  const double a[n_wts] = {3.9355079629947969884346551941969960e-01L,
                                           4.7979065808897448654107733982929214e-01L,
                                           5.1003445645828061436081405648347852e-03L,
                                           2.6597620190330158952732822450744488e-01L,
                                           2.8536418538696461608233522814483715e-01L,
                                           1.3723536747817085036455583801851025e-01L};

                  const double b[n_wts] = {0.,
                                           0.,
                                           5.6817155788572446538150614865768991e-02L,
                                           1.2539956353662088473247489775203396e-01L,
                                           1.2409970153698532116262152247041742e-02L,
                                           5.2792057988217708934207928630851643e-02L};

                  const unsigned int permutation_ids[n_wts] = {3, 3, 6, 6, 6, 6}; // 30 total points

                  dunavant_rule2(wts, a, b, permutation_ids, n_wts);

                  return;
                }
            } // end of switch(order)
        }     // end of case(quad)

        //---------------------------------------------
        // Unsupported type
      default: assert("Element type not supported! " && 0);
    }
}


/* -----------------------------------------------------------------------
 *   3D rules
 * ----------------------------------------------------------------------- */

void
QGenerator::conical_rule_tet()
{
  // Be sure the underlying rule object was built with the same dimension as the
  // rule we are about to construct.
  assert(cell_dim(_type) == 3);

  QGenerator gauss1D;
  gauss1D.init(CELL_TYPE_EDGE, _order);
  gauss1D.init_1d_gauss();

  QGenerator jacA1D;
  jacA1D.init(CELL_TYPE_EDGE, _order);
  jacA1D.init_1d_jacobi(1);

  QGenerator jacB1D;
  jacB1D.init(CELL_TYPE_EDGE, _order);
  jacB1D.init_1d_jacobi(2);

  // The Gauss rule needs to be scaled to [0,1]
  double new_range[] = {0, 1};
  gauss1D.scale(new_range);

  // Now construct the points and weights for the conical product rule.

  // All rules should have the same number of points
  assert(gauss1D.n_points() == jacA1D.n_points());
  assert(jacA1D.n_points() == jacB1D.n_points());

  // Save the number of points as a convenient variable
  const unsigned int np = gauss1D.n_points();

  // All rules should be between x=0 and x=1
  assert(gauss1D.qp(0)(0) > 0.0);
  assert(gauss1D.qp(np - 1)(0) < 1.0);
  assert(jacA1D.qp(0)(0) > 0.0);
  assert(jacA1D.qp(np - 1)(0) < 1.0);
  assert(jacB1D.qp(0)(0) > 0.0);
  assert(jacB1D.qp(np - 1)(0) < 1.0);

  // Resize the points and weights vectors
  _points.resize(np * np * np);
  _weights.resize(np * np * np);

  // Compute the conical product
  unsigned int gp = 0;
  for(unsigned int i = 0; i < np; i++)
    for(unsigned int j = 0; j < np; j++)
      for(unsigned int k = 0; k < np; k++)
        {
          // t[k];
          _points[gp](0) = jacB1D.qp(k)(0);
          // s[j]*(1.-t[k]);
          _points[gp](1) = jacA1D.qp(j)(0) * (1. - jacB1D.qp(k)(0));
          // r[i]*(1.-s[j])*(1.-t[k]);
          _points[gp](2) = gauss1D.qp(i)(0) * (1. - jacA1D.qp(j)(0)) * (1. - jacB1D.qp(k)(0));
          // A[i]*B[j]*C[k];
          _weights[gp] = gauss1D.w(i) * jacA1D.w(j) * jacB1D.w(k);
          gp++;
        }
}

void
QGenerator::tensor_product_hex()
{
  QGenerator q1D;
  q1D.init(CELL_TYPE_EDGE, _order);
  q1D.init_1d_gauss();

  const unsigned int np = q1D.n_points();

  _points.resize(np * np * np);

  _weights.resize(np * np * np);

  unsigned int q = 0;

  for(unsigned int k = 0; k < np; k++)
    for(unsigned int j = 0; j < np; j++)
      for(unsigned int i = 0; i < np; i++)
        {
          _points[q](0) = q1D.qp(i)(0);
          _points[q](1) = q1D.qp(j)(0);
          _points[q](2) = q1D.qp(k)(0);

          _weights[q] = q1D.w(i) * q1D.w(j) * q1D.w(k);

          q++;
        }
}

void
QGenerator::tensor_product_prism()
{
  QGenerator q1D;
  q1D.init(CELL_TYPE_EDGE, _order);
  q1D.init_1d_gauss();

  QGenerator q2D;
  q2D.init(CELL_TYPE_TRI, _order);
  q2D.init_2d_gauss();

  const unsigned int n_points1D = q1D.n_points();
  const unsigned int n_points2D = q2D.n_points();

  _points.resize(n_points1D * n_points2D);
  _weights.resize(n_points1D * n_points2D);

  unsigned int q = 0;

  for(unsigned int j = 0; j < n_points1D; j++)
    for(unsigned int i = 0; i < n_points2D; i++)
      {
        _points[q](0) = q2D.qp(i)(0);
        _points[q](1) = q2D.qp(i)(1);
        _points[q](2) = q1D.qp(j)(0);

        _weights[q] = q2D.w(i) * q1D.w(j);

        q++;
      }
}

// Builds and scales a Gauss rule and a Jacobi rule.
// Then combines them to compute points and weights
// of a 3D conical product rule for the Pyramid.  The integral
// over the reference Tet can be written (in LaTeX notation) as:
//
// If := \int_0^1 dz \int_{-(1-z)}^{(1-z)} dy \int_{-(1-z)}^{(1-z)} f(x,y,z) dx (1)
//
// (Imagine a stack of infinitely thin squares which decrease in size as
//  you approach the apex.)  Under the transformation of variables:
//
// z=w
// y=(1-z)v
// x=(1-z)u,
//
// The Jacobian determinant of this transformation is |J|=(1-w)^2, and
// the integral itself is transformed to:
//
// If = \int_0^1 (1-w)^2 dw \int_{-1}^{1} dv \int_{-1}^{1} f((1-w)u, (1-w)v, w) du (2)
//
// The integral can now be approximated by the product of three 1D quadrature rules:
// A Jacobi rule with alpha==2, beta==0 in w, and Gauss rules in v and u.  In this way
// we can obtain 3D rules to any order for which the 1D rules exist.
void
QGenerator::conical_rule_pyramid()
{
  // Be sure the underlying rule object was built with the same dimension as the
  // rule we are about to construct.
  assert(cell_dim(_type) == 3);

  QGenerator gauss1D;
  gauss1D.init(CELL_TYPE_EDGE, _order);
  gauss1D.init_1d_gauss();

  QGenerator jac1D;
  jac1D.init(CELL_TYPE_EDGE, _order);
  jac1D.init_1d_jacobi(2);

  // These rules should have the same number of points
  assert(gauss1D.n_points() == jac1D.n_points());

  // Save the number of points as a convenient variable
  const unsigned int np = gauss1D.n_points();

  // Resize the points and weights vectors
  _points.resize(np * np * np);
  _weights.resize(np * np * np);

  // Compute the conical product
  unsigned int q = 0;
  for(unsigned int i = 0; i < np; ++i)
    for(unsigned int j = 0; j < np; ++j)
      for(unsigned int k = 0; k < np; ++k, ++q)
        {
          const double xi = gauss1D.qp(i)(0);
          const double yj = gauss1D.qp(j)(0);
          const double zk = jac1D.qp(k)(0);

          _points[q](0) = (1. - zk) * xi;
          _points[q](1) = (1. - zk) * yj;
          _points[q](2) = zk;
          _weights[q] = gauss1D.w(i) * gauss1D.w(j) * jac1D.w(k);
        }
}


void
QGenerator::init_3d_gauss()
{

  switch(_type)
    {

        //---------------------------------------------
        // Hex quadrature rules
      case CELL_TYPE_HEX:
        {
          // We compute the 3D quadrature rule as a tensor
          // product of the 1D quadrature rule.
          tensor_product_hex();
          return;
        }


        //---------------------------------------------
        // Tetrahedral quadrature rules
      case CELL_TYPE_TET:
        {
          switch(_order)
            {
                // Taken from pg. 222 of "The finite element method," vol. 1
                // ed. 5 by Zienkiewicz & Taylor
              case 0:
              case 1:
                {
                  // Exact for linears
                  _points.resize(1);
                  _weights.resize(1);


                  _points[0](0) = .25;
                  _points[0](1) = .25;
                  _points[0](2) = .25;

                  _weights[0] = double(1) / 6;

                  return;
                }
              case 2:
                {
                  // Exact for quadratics
                  _points.resize(4);
                  _weights.resize(4);


                  const double a = .585410196624969;
                  const double b = .138196601125011;

                  _points[0](0) = a;
                  _points[0](1) = b;
                  _points[0](2) = b;

                  _points[1](0) = b;
                  _points[1](1) = a;
                  _points[1](2) = b;

                  _points[2](0) = b;
                  _points[2](1) = b;
                  _points[2](2) = a;

                  _points[3](0) = b;
                  _points[3](1) = b;
                  _points[3](2) = b;



                  _weights[0] = double(1) / 24;
                  _weights[1] = _weights[0];
                  _weights[2] = _weights[0];
                  _weights[3] = _weights[0];

                  return;
                }


              // Can be found in the class notes
              // http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps
              // by Flaherty.
              //
              // Caution: this rule has a negative weight and may be
              // unsuitable for some problems.
              // Exact for cubics.
              //
              // Note: Keast (see ref. elsewhere in this file) also gives
              // a third-order rule with positive weights, but it contains points
              // on the ref. elt. boundary, making it less suitable for FEM calculations.
              case 3:
                {
                  // If a rule with postive weights is required, the 2x2x2 conical
                  // product rule is third-order accurate and has less points than
                  // the next-available positive-weight rule at FIFTH order.
                  conical_rule_tet();
                  return;
                }


              // Originally a Keast rule,
              //    Patrick Keast,
              //    Moderate Degree Tetrahedral QGenerator Formulas,
              //    Computer Methods in Applied Mechanics and Engineering,
              //    Volume 55, Number 3, May 1986, pages 339-348.
              //
              // Can also be found the class notes
              // http://www.cs.rpi.edu/~flaherje/FEM/fem6.ps
              // by Flaherty.
              //
              // Caution: this rule has a negative weight and may be
              // unsuitable for some problems.
              case 4:


              // Walkington's fifth-order 14-point rule from
              // "QGenerator on Simplices of Arbitrary Dimension"
              //
              // We originally had a Keast rule here, but this rule had
              // more points than an equivalent rule by Walkington and
              // also contained points on the boundary of the ref. elt,
              // making it less suitable for FEM calculations.
              case 5:
                {
                  _points.resize(14);
                  _weights.resize(14);

                  // permutations of these points and suitably-modified versions of
                  // these points are the quadrature point locations
                  const double a[3] = {0.31088591926330060980,   // a1 from the paper
                                       0.092735250310891226402,  // a2 from the paper
                                       0.045503704125649649492}; // a3 from the paper

                  // weights.  a[] and wt[] are the only floating-point inputs required
                  // for this rule.
                  const double wt[3] = {0.018781320953002641800,   // w1 from the paper
                                        0.012248840519393658257,   // w2 from the paper
                                        0.0070910034628469110730}; // w3 from the paper

                  // The first two sets of 4 points are formed in a similar manner
                  for(unsigned int i = 0; i < 2; ++i)
                    {
                      // Where we will insert values into _points and _weights
                      const unsigned int offset = 4 * i;

                      // Stuff points and weights values into their arrays
                      const double b = 1. - 3. * a[i];

                      // Here are the permutations.  Order of these is not important,
                      // all have the same weight
                      _points[offset + 0] = Point(a[i], a[i], a[i]);
                      _points[offset + 1] = Point(a[i], b, a[i]);
                      _points[offset + 2] = Point(b, a[i], a[i]);
                      _points[offset + 3] = Point(a[i], a[i], b);

                      // These 4 points all have the same weights
                      for(unsigned int j = 0; j < 4; ++j) _weights[offset + j] = wt[i];
                    } // end for


                  {
                    // The third set contains 6 points and is formed a little differently
                    const unsigned int offset = 8;
                    const double b = 0.5 * (1. - 2. * a[2]);

                    // Here are the permutations.  Order of these is not important,
                    // all have the same weight
                    _points[offset + 0] = Point(b, b, a[2]);
                    _points[offset + 1] = Point(b, a[2], a[2]);
                    _points[offset + 2] = Point(a[2], a[2], b);
                    _points[offset + 3] = Point(a[2], b, a[2]);
                    _points[offset + 4] = Point(b, a[2], b);
                    _points[offset + 5] = Point(a[2], b, b);

                    // These 6 points all have the same weights
                    for(unsigned int j = 0; j < 6; ++j) _weights[offset + j] = wt[2];
                  }

                  // for(uint i = 0 ; i < 14 ; i++)
                  // 	printf("%lf \n", _weights[i]);
                  // assert(0);
                  return;
                }

              // This rule is originally from Keast:
              //    Patrick Keast,
              //    Moderate Degree Tetrahedral QGenerator Formulas,
              //    Computer Methods in Applied Mechanics and Engineering,
              //    Volume 55, Number 3, May 1986, pages 339-348.
              //
              // It is accurate on 6th-degree polynomials and has 24 points
              // vs. 64 for the comparable conical product rule.
              //
              // Values copied 24th June 2008 from:
              // http://people.scs.fsu.edu/~burkardt/f_src/keast/keast.f90
              case 6:

              // Keast's 31 point, 7th-order rule contains points on the reference
              // element boundary, so we've decided not to include it here.
              //
              // Keast's 8th-order rule has 45 points.  and a negative
              // weight, so if you've explicitly disallowed such rules
              // you will fall through to the conical product rules
              // below.
              case 7:
              case 8:

              // Fall back on Grundmann-Moller or Conical Product rules at high orders.
              default:
                {
                  // The following quadrature rules are generated as
                  // conical products.  These tend to be non-optimal
                  // (use too many points, cluster points in certain
                  // regions of the domain) but they are quite easy to
                  // automatically generate using a 1D Gauss rule on
                  // [0,1] and two 1D Jacobi-Gauss rules on [0,1].
                  conical_rule_tet();
                  return;
                }
            } // end of switch(order)

        } // end case CELL_TYPE_TET


        //---------------------------------------------
        // CELL_TYPE_WEDGE quadrature rules
      case CELL_TYPE_WEDGE:
        {
          // We compute the 3D quadrature rule as a tensor
          // product of the 1D quadrature rule and a 2D
          // triangle quadrature rule
          tensor_product_prism();

          return;
        }


        //---------------------------------------------
        // Pyramid
      case CELL_TYPE_PYRAMID:
        {
          conical_rule_pyramid();
          return;
        }



        //---------------------------------------------
        // Unsupported type
      default: assert("ERROR: Unsupported 3D cell type: " && 0);
    }
}