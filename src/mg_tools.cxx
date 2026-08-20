#include "mg_tools.hxx"
#include "grid_iterators.hxx"

namespace mg
{

Context::Context()
: parameters()
, grids()
, transfer_infos()
, AA()
, smoothers()
, uu()
, rr()
, du()
, _is_set_smoothers_called(false)
, _is_set_grids_called(false)
, _is_set_parameters_called(false)
, _is_set_up_transfer_info_called(false)
, _is_set_up_discrete_equations_called(false)
, _tolerance_factor(1e-10)
, _max_iterations(20)
{
}

void
Context::verify() const
{

  assert(n_levels() >= 1);
  assert((int)grids.size() == n_levels());
  assert((int)transfer_infos.size() == n_levels() - 1);
  assert((int)AA.size() == n_levels());
  assert((int)rr.size() == n_levels());
  assert((int)uu.size() == n_levels());
  assert((int)du.size() == n_levels());
  assert((int)ff.size() == n_levels());
  //
  assert(_is_set_grids_called);
  assert(_is_set_parameters_called);
  assert(_is_set_up_transfer_info_called);
  assert(_is_set_up_discrete_equations_called);
  assert(_is_set_smoothers_called);
}

void
Context::set_parameters(const Parameters & parameters_in)
{
  parameters = parameters_in;
  _is_set_parameters_called = true;
  //
  grids.resize(n_levels());
  transfer_infos.resize(n_levels() - 1);
  //
  AA.resize(n_levels());
  smoothers.resize(n_levels());
  uu.resize(n_levels());
  rr.resize(n_levels());
  du.resize(n_levels());
  ff.resize(n_levels());
}

void
Context::set_grids(std::vector<Grid_connectivity_context> & grids_in)
{
  assert(_is_set_parameters_called);
  assert((int)grids_in.size() == n_levels());
  //
  for(int i = 0; i < n_levels(); ++i)
    {
      grids[i].swap(grids_in[i]);
    }
  //
  for(int i = 1; i < n_levels(); ++i)
    {
      assert(grids[i - 1].n_vertices() > grids[i].n_vertices());
    }

  _is_set_grids_called = true;
}

void
Context::set_up_transfer_info()
{
  int ierr;
  assert(_is_set_grids_called);
  for(int level = 0; level < n_levels() - 1; ++level)
    {
      // tr[i] locates i on i+1
      Point_locator locator(grids[level + 1]);
      ierr = locator.bfs_mesh_search(grids[level], transfer_infos[level]);
      assert(ierr == Point_locator::LOCATOR_SUCCESS);
    }
  _is_set_up_transfer_info_called = true;
}

void
Context::set_up_discrete_equations(const Poisson_assembler & assembler)
{
  assert(_is_set_grids_called);

  for(int level = 0; level < n_levels(); ++level)
    {
      const int n_dofs = assembler.get_n_dof(grids[level]);
      assembler.get_nonzero_structure(grids[level], AA[level]);
      for(int b = 0; b < n_buffers; ++b)
        {
          uu[level][b].init(n_dofs);
          ff[level][b].init(n_dofs);
          rr[level][b].init(n_dofs);
          du[level][b].init(n_dofs);
        }
      // The higher level ff's are actually ignored
      // only level 0 matters
      assembler.assemble_system(grids[level], AA[level], ff[level][ZERO]);
    }
  _is_set_up_discrete_equations_called = true;
}


void
Context::set_smoothers(std::vector<Preconditioner *> & precons_in)
{
  assert(_is_set_up_discrete_equations_called);
  for(int i = 0; i < n_levels(); ++i)
    {
      smoothers[i].reset(precons_in[i]);
      precons_in[i] = NULL;
      smoothers[i]->form(AA[i]);
    }
  _is_set_smoothers_called = true;
}

void
Context::set_smoothers(const Preconditioner & precon_in)
{
  assert(_is_set_up_discrete_equations_called);
  for(int i = 0; i < n_levels(); ++i)
    {
      smoothers[i].reset(precon_in.clone());
      smoothers[i]->form(AA[i]);
    }
  _is_set_smoothers_called = true;
}

void
Context::prolongate_dofs(const int level_from, const int level_to, const Numeric_array & in, Numeric_array & out)
{
  assert(level_to == level_from - 1);
  assert(level_to >= 0);
  assert(level_from < n_levels());
  verify();

  Context::interpolate_impl(transfer_infos[level_to], grids[level_from], in, grids[level_to], out);
}

void
Context::restrict_dofs(const int level_from, const int level_to, const Numeric_array & in, Numeric_array & out)
{
  assert(level_to == level_from + 1);
  assert(level_from >= 0);
  assert(level_to < n_levels());
  verify();

  Context::interpolate_transpose_impl(transfer_infos[level_from], grids[level_from], in, grids[level_to], out);
}


void
Context::interpolate_impl(const Transfer_info & tr,
                          const Grid_connectivity_context & gin,
                          const Numeric_array & in,
                          const Grid_connectivity_context & gout,
                          Numeric_array & out)
{
  Grid_iterators giterin(gin);
  // TODO: make these faster by using an actuall Prolongation sparse matrix

  // Prolongation is interpolation
  // so tr should know where out is located on in.
  assert((int)in.size() == gin.n_vertices());
  assert((int)out.size() == gout.n_vertices());
  assert(tr.owner_id.size() == out.size());

  for(int i = 0; i < gout.n_vertices(); ++i)
    {
      const Grid_iterators::Triangle triowner = giterin.triangle_iterator(tr.owner_id[i]);

      if(triowner.is_dummy())
        {
          const int vownerid0 = triowner.vertex_from_offset(1).id();
          const int vownerid1 = triowner.vertex_from_offset(2).id();
          const double xi0 = 1 - tr.xi[2 * i + 0]; // TODO this has to be fixed in the API (it should be xi not 1-xi)
          const double xi1 = tr.xi[2 * i + 0];
          assert(xi0 > -1e-10);
          assert(xi1 > -1e-10);
          assert(xi0 < 1 + 1e-10);
          assert(xi1 < 1 + 1e-10);
          assert(tr.xi[2 * i + 1] < -100);
          out[i] = xi0 * in[vownerid0] + xi1 * in[vownerid1];
        }
      else
        {
          const int vownerid0 = triowner.vertex_from_offset(0).id();
          const int vownerid1 = triowner.vertex_from_offset(1).id();
          const int vownerid2 = triowner.vertex_from_offset(2).id();
          const double xi0 = tr.xi[2 * i + 0];
          const double xi1 = tr.xi[2 * i + 1];
          const double xi2 = 1 - xi0 - xi1;
          assert(xi0 > -1e-10);
          assert(xi1 > -1e-10);
          assert(xi2 > -1e-10);
          assert(xi0 < 1 + 1e-10);
          assert(xi1 < 1 + 1e-10);
          assert(xi2 < 1 + 1e-10);
          out[i] = xi0 * in[vownerid0] + xi1 * in[vownerid1] + xi2 * in[vownerid2];
        }
    }
} // All done

void
Context::interpolate_transpose_impl(const Transfer_info & tr,
                                    const Grid_connectivity_context & gin,
                                    const Numeric_array & in,
                                    const Grid_connectivity_context & gout,
                                    Numeric_array & out)
{
  Grid_iterators giterout(gout);

  // TODO: make these faster by using an actuall Prolongation sparse matrix
  //  This needs set zero
  out.set_zero();

  // Numeric_array interim (out.size());
  // interim.set_zero();

  // Restriction is the transpose of  interpolation
  // tr should know where in is located on out.
  assert((int)in.size() == gin.n_vertices());
  assert((int)out.size() == gout.n_vertices());
  assert(tr.owner_id.size() == in.size());

  for(int i = 0; i < gin.n_vertices(); ++i)
    {
      const Grid_iterators::Triangle triowner = giterout.triangle_iterator(tr.owner_id[i]);

      if(triowner.is_dummy())
        {
          const int vownerid0 = triowner.vertex_from_offset(1).id();
          const int vownerid1 = triowner.vertex_from_offset(2).id();
          const double xi0 = 1 - tr.xi[2 * i + 0]; // TODO this has to be fixed in the API (it should be xi not 1-xi)
          const double xi1 = tr.xi[2 * i + 0];
          assert(xi0 > -1e-10);
          assert(xi1 > -1e-10);
          assert(xi0 < 1 + 1e-10);
          assert(xi1 < 1 + 1e-10);
          assert(tr.xi[2 * i + 1] < -100);
          out[vownerid0] += xi0 * in[i];
          out[vownerid1] += xi1 * in[i];
          //
          // interim[vownerid0] += xi0;
          // interim[vownerid1] += xi1;
        }
      else
        {
          const int vownerid0 = triowner.vertex_from_offset(0).id();
          const int vownerid1 = triowner.vertex_from_offset(1).id();
          const int vownerid2 = triowner.vertex_from_offset(2).id();
          const double xi0 = tr.xi[2 * i + 0];
          const double xi1 = tr.xi[2 * i + 1];
          const double xi2 = 1 - xi0 - xi1;
          assert(xi0 > -1e-10);
          assert(xi1 > -1e-10);
          assert(xi2 > -1e-10);
          assert(xi0 < 1 + 1e-10);
          assert(xi1 < 1 + 1e-10);
          assert(xi2 < 1 + 1e-10);
          out[vownerid0] += xi0 * in[i];
          out[vownerid1] += xi1 * in[i];
          out[vownerid2] += xi2 * in[i];
          //
          // interim[vownerid0] += xi0;
          // interim[vownerid1] += xi1;
          // interim[vownerid2] += xi2;
        }
    }

  // for (int i = 0 ; i < (int) out.size() ; ++i)
  //{
  //  out[i] /= interim[i];
  //}

} // All done

double
Context::smooth(const int level, const int b, const int n_times)
{
  double rnorm;
  verify();

  // u = u + P(b-Au):
  for(int i = 0; i < n_times; ++i)
    {
      // r = b - Au
      rnorm = eval_residual(AA[level], uu[level][b], ff[level][b], rr[level][b]);

      // du <-  P(r)
      smoothers[level]->apply(rr[level][b], du[level][b]);

      // u = u + du
      uu[level][b].add(du[level][b]);
    }

  return rnorm;
}

double
Context::eval_residual(const Sparse_matrix & A, const Numeric_array & u, const Numeric_array & f, Numeric_array & r)
{
  r.set_equal(f);
  bridson::multiply_scale_and_add(A, u, -1, r);
  return r.norm_l2();
}


double
Context::cycle(const int level, const int b)
{
  int tau;
  double rnorm;

  verify();

  switch(parameters.cycle_type)
    {
      case Parameters::CYCLE_TYPE_V: tau = 1; break;
      case Parameters::CYCLE_TYPE_W: tau = 2; break;
      default: tau = -1; assert(0);
    }

  if(level == n_levels() - 1)
    {
      rnorm = this->smooth(level, b, n_smoothing(level));
      // printf(" Level %d coarse-solve: %g \n", level, rnorm);
    }
  else
    {
      assert(level >= 0);
      assert(level < n_levels() - 1);

      // Perform pre-smoothing
      rnorm = this->smooth(level, b, n_smoothing(level));
      // printf(" Level %d pre-smooth: %g \n", level, rnorm);

      // Find the residual
      rnorm = eval_residual(AA[level], uu[level][b], ff[level][b], rr[level][b]);

      // Move the residual to the higher level
      this->restrict_dofs(level, level + 1, rr[level][b], ff[level + 1][b]);

      // Set the init uu on higher level to zero.
      uu[level + 1][b].set_zero();

      // Perform cycles on the solution on the upper level
      for(int i = 0; i < tau; ++i)
        {
          rnorm = this->cycle(level + 1, b);
          // printf(" Level %d cycle: %g \n", level, rnorm);
        }

      // Move the correction from the upper level to the current level
      this->prolongate_dofs(level + 1, level, uu[level + 1][b], du[level][b]);

      // and add it to the solution on currecnt level
      uu[level][b].add_scaled(1, du[level][b]);

      // Post smoothing
      rnorm = this->smooth(level, b, n_smoothing(level));
      // printf(" Level %d post-smooth: %g \n", level, rnorm);
    }

  // if(level == 0) printf("After cycle: %g \n", rnorm);
  return rnorm;

} // All done



} // End of namespace mg