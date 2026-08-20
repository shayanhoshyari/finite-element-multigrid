#include "poisson_equation.hxx"
#include "fe.hxx"
#include "grid_iterators.hxx"


int
Poisson_assembler::get_n_dof(const Grid_connectivity_context & grid) const
{
  return grid.n_vertices();
}

void
Poisson_assembler::get_nonzero_structure(const Grid_connectivity_context & grid,
                                         std::vector<int> & xadj,
                                         std::vector<int> & adj) const
{
  Grid_iterators giter(grid);

  //
  // Count nnz per row
  //
  xadj.resize(0);
  xadj.reserve(giter.n_vertices());
  xadj.push_back(0);
  for(int i = 0; i < giter.n_vertices(); ++i)
    {
      Grid_iterators::Vertex_umbrella umb = giter.vertex_umbrella_iterator(i);
      xadj.push_back(xadj.back());
      ++xadj.back(); // vi
      do
        {
          ++xadj.back(); // vj
        }
      while(umb.advance());
    }

  // Now records
  adj.resize(0);
  adj.reserve(xadj.back());
  for(int i = 0; i < giter.n_vertices(); ++i)
    {
      Grid_iterators::Vertex_umbrella umb = giter.vertex_umbrella_iterator(i);
      Grid_iterators::Vertex vertex_i = giter.vertex_iterator(i);

      adj.push_back(vertex_i.id());
      do
        {
          Grid_iterators::Vertex vertex_j = umb.half_edge().origin();
          adj.push_back(vertex_j.id());
        }
      while(umb.advance());
    }
}

void
Poisson_assembler::get_nonzero_structure(const Grid_connectivity_context & grid, bridson::FixedSparseMatrix & lhs) const
{
  std::vector<int> adj, xadj;
  this->get_nonzero_structure(grid, xadj, adj);
  lhs.set_adjacency(xadj, adj);
}

void
Poisson_assembler::assemble_system(const Grid_connectivity_context & grid,
                                   bridson::FixedSparseMatrix & lhs,
                                   Numeric_array & rhs) const
{
  const int n_dof = get_n_dof(grid);
  const Grid_iterators giter(grid);

  //
  // Ensure that the sizes are correct and set both LHS and RHS to zero
  //
  assert(lhs.n == n_dof);
  lhs.set_zero();
  if(rhs.size())
    {
      assert((int)rhs.size() == n_dof);
      rhs.set_zero();
    }

  //
  // Get access to the quadrature information
  // For now, set qorder to 3. In reality, 1 is enough for LHS but RHS needs higher.
  //
  const int qorder = 3;

  // Triangles
  FE_triangle fe_tri(qorder);
  const std::vector<bridson::Vec2d> & qp_tri = fe_tri.qp();
  const std::vector<double> & JxW_tri = fe_tri.JxW();
  const std::vector<std::vector<double>> & phi_tri = fe_tri.phi();
  const std::vector<std::vector<bridson::Vec2d>> & dphi_tri = fe_tri.dphi();

  // segments
  FE_segment fe_seg(qorder);
  const std::vector<bridson::Vec2d> & qp_seg = fe_seg.qp();
  const std::vector<double> & JxW_seg = fe_seg.JxW();
  const std::vector<std::vector<double>> & phi_seg = fe_seg.phi();
  // const std::vector<std::vector<bridson::Vec2d>> & dphi_seg = fe_seg.dphi();

  const bridson::Vec4d stiffness = problem_data().stiffness();

  //
  // Loop over all internal
  //
  for(int triid = 0; triid < giter.n_real_triangles(); ++triid)
    {
      bridson::Vec3d fe(0);
      bridson::Vec<9, double> Ke(0);
      bridson::Vec3i dof_ids(0);

      const Grid_iterators::Triangle triiter = giter.triangle_iterator(triid);
      const Grid_iterators::Vertex viter0 = triiter.vertex_from_offset(0);
      const Grid_iterators::Vertex viter1 = triiter.vertex_from_offset(1);
      const Grid_iterators::Vertex viter2 = triiter.vertex_from_offset(2);

      dof_ids[0] = viter0.id();
      dof_ids[1] = viter1.id();
      dof_ids[2] = viter2.id();

      const bridson::Vec2d v0pos(viter0.x(), viter0.y());
      const bridson::Vec2d v1pos(viter1.x(), viter1.y());
      const bridson::Vec2d v2pos(viter2.x(), viter2.y());

      fe_tri.reinit(v0pos, v1pos, v2pos);

      const unsigned n_phi = phi_tri.size(); // fancy assness, this is obviously 3.

      for(unsigned int qp_id = 0; qp_id < qp_tri.size(); ++qp_id)
        {
          const double ff = problem_data().right_hand_side(qp_tri[qp_id]);
          for(unsigned i = 0; i < n_phi; ++i)
            {
              for(unsigned j = 0; j < n_phi; ++j)
                {
                  // TODO: incase of non-symmetric K check that this is right.
                  Ke[i * n_phi + j] +=
                    JxW_tri[qp_id] * bridson::dot(dphi_tri[i][qp_id], bridson::mat_prod(stiffness, dphi_tri[j][qp_id]));
                }
              fe[i] += JxW_tri[qp_id] * ff * phi_tri[i][qp_id];
            }
        } //  End of quadrature


      // Contribute to the global system
      if(rhs.size()) rhs.add_values(dof_ids, fe);
      lhs.add_values(dof_ids, dof_ids, Ke);
    } // All done assembling internal elements

  //
  // Loop over all boundary elements
  // (Or it can be done at same time with the internal elements).
  //
  for(int bedgeid = 0; bedgeid < giter.n_boundary_edges(); ++bedgeid)
    {
      bridson::Vec2d fe(0);
      bridson::Vec<4, double> Ke(0);
      bridson::Vec2i dof_ids(0);

      const Grid_iterators::Half_edge heiter = giter.boundary_half_edge_iterator(bedgeid);
      const Grid_iterators::Vertex viter0 = heiter.origin();
      const Grid_iterators::Vertex viter1 = heiter.twin().origin();
      const int btag = heiter.triangle().tag();

      dof_ids[0] = viter0.id();
      dof_ids[1] = viter1.id();

      const bridson::Vec2d v0pos(viter0.x(), viter0.y());
      const bridson::Vec2d v1pos(viter1.x(), viter1.y());

      fe_seg.reinit(v0pos, v1pos);

      const unsigned n_phi = phi_seg.size(); // fancy assness, this is obviously 2.

      for(unsigned int qp_id = 0; qp_id < qp_seg.size(); ++qp_id)
        {
          double qq, tt;
          problem_data().boundary_condition_terms(qp_seg[qp_id], btag, qq, tt);
          for(unsigned i = 0; i < n_phi; ++i)
            {
              for(unsigned j = 0; j < n_phi; ++j)
                {
                  // TODO: incase of non-symmetric K check that this is right.
                  Ke[i * n_phi + j] += JxW_seg[qp_id] * qq * phi_seg[i][qp_id] * phi_seg[j][qp_id];
                }
              fe[i] += JxW_seg[qp_id] * tt * phi_tri[i][qp_id];
            }
        } //  End of quadrature


      // Contribute to the global system
      if(rhs.size()) rhs.add_values(dof_ids, fe);
      lhs.add_values(dof_ids, dof_ids, Ke);
    } // All done assembling boundary elements
}

void
Poisson_assembler::find_exact_dofs(const Grid_connectivity_context & grid, Numeric_array & dofs) const
{
  assert(problem_data().has_exact_solution());
  assert((int)dofs.size() == get_n_dof(grid));

  const Grid_iterators giter(grid);


  for(int vid = 0; vid < giter.n_vertices(); ++vid)
    {
      Grid_iterators::Vertex viter = giter.vertex_iterator(vid);
      const bridson::Vec2d pos(viter.x(), viter.y());
      dofs[vid] = problem_data().exact_solution(pos);
    }
}
