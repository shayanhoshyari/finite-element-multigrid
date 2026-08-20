/*
 * grid_io.hxx
 *
 * A variety of functionality for reading and writing meshes in different
 * formats.
 */

#ifndef GRID_IO_HXX_IS_INCLUDED
#define GRID_IO_HXX_IS_INCLUDED

#include <string>

#include "grid_connectivity_context.hxx"

class Grid_IO
{
  // pointer to the mesh that we are working on.
  Grid_connectivity_context & _grid;
  bool _have_written_cell_data_header;
  bool _have_written_vert_data_header;
  int _writing_mode;

public:
  enum
  {
    MODE_CELL = 0,
    MODE_VERTEX = 1
  };

  // Trivial constructor
  Grid_IO(Grid_connectivity_context & mesh_in);
  Grid_connectivity_context & grid() { return _grid; }
  const Grid_connectivity_context & grid() const { return _grid; }

  // Functions to read a mesh.
  void read_auto(const std::string &);
  void read_obj(const std::string &);
  void read_off(const std::string &);
  void read_gmsh(const std::string &);
  void read_vtk(const std::string &);
  void read_triangle(const std::string &);

  // Functions to write a mesh
  void write_off(const std::string &) const;
  void write_vtk(const std::string &);
  void write_vtk(FILE *);
  void write_vtk_vert_header(FILE *);
  void write_vtk_cell_header(FILE *);
  void write_vtk_data(FILE *, std::vector<int> &, const std::string name, const bool is_vector=false );
  void write_vtk_data(FILE *, std::vector<double> &, const std::string name, const bool is_vector=false);
  void print_data_info(FILE * fl = stdout) const;
  void print_iter_info(FILE * fl = stdout) const;
  // void write_vtk_dual(const std::string & fname);

  // Open a file
  static FILE * open_file(const std::string & fname, const std::string & format, const std::string & mode);
};

#endif /* GRID_IO_HXX_IS_INCLUDED */