#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>

#include "grid_io.hxx"
#include "grid_iterators.hxx"


using std::string;
using std::vector;

// clang-format off
#define io_ensure(cmd, msg) do{const bool success=(cmd); assert(success &&msg);  if(!success){printf("[ERROR] %s", msg); throw;} }while(0)
// clang-format on

FILE *
Grid_IO::open_file(const string & fname, const string & format, const string & mode)
{
  // open the file
  FILE * fl;
  string final_fname;

  if(fname.find(format) == fname.length() - format.length())
    {
      final_fname = fname;
    }
  else
    {
      final_fname = fname + format;
    }

  fl = fopen(final_fname.c_str(), mode.c_str());
  if(!fl)
    {
      std::cout << "[Error] Could not find the file " << final_fname << "." << std::endl;
      assert(0);
      throw "Bad input file";
    }

  return fl;
}

Grid_IO::Grid_IO(Grid_connectivity_context & mesh_in)
: _grid(mesh_in)
, _have_written_cell_data_header(false)
, _have_written_vert_data_header(false)
, _writing_mode(-1)
{
}

void
Grid_IO::read_auto(const string & fname)
{
  auto has_extension = [](const string & f, const string & e) -> bool { return f.find(e) == f.length() - e.length(); };

  if(has_extension(fname, ".obj"))
    {
      this->read_obj(fname);
    }
  else if(has_extension(fname, ".off"))
    {
      this->read_off(fname);
    }
  else if(has_extension(fname, ".msh"))
    {
      this->read_gmsh(fname);
    }
  else
    {
      std::cout << "Unrecognized format " << fname << std::endl;
      throw "Unsupported input file format";
    }
}

void
Grid_IO::read_obj(const string & fname)
{
  // open the file
  FILE * fl = open_file(fname, ".obj", "r");

  /*
   * read the file
   */
  vector<double> vcoord;
  vector<int> ftovert;

  // Read each line
  int face_data_mode = -1;
  for(;;)
    {

      int c;
      c = fgetc(fl);

      if(c == EOF)
        {
          break;
        }
      else if(c == (int)'v')
        {
          c = fgetc(fl);
          if(c == (int)' ')
            {
              double x, y, z;
              const int success = fscanf(fl, "%lf %lf %lf", &x, &y, &z);
              assert(success == 3);
              vcoord.push_back(x);
              vcoord.push_back(y);
              vcoord.push_back(z);
            }
        }
      else if(c == (int)'f')
        {

          // read the face data in a single line
          char line[1000];
          char * success = fgets(line, 1000, fl);
          assert(success == line);

          // find how the verts are written in a line
          // if it is not already found.
          if(face_data_mode < 0)
            {
              // read the first word in the string
              char word[100];
              sscanf(line, "%s", word);
              // printf("word: %s \n", word);
              const int n_char = strlen(word);
              int n_slash = 0;

              // count the number of slashes in there
              // and decide the format based on that.
              for(int i = 0; i < n_char; i++)
                {
                  if(word[i] == '/')
                    {
                      if(word[i + 1] == '/')
                        {
                          face_data_mode = 3;
                          break;
                        }
                      else
                        {
                          n_slash++;
                        }
                    }
                }

              if(n_slash == 2)
                face_data_mode = 4;
              else if(n_slash == 1)
                face_data_mode = 2;
              else if(face_data_mode != 3)
                face_data_mode = 1;

              std::cerr << "MeshIO::read_obj: mode= " << face_data_mode << std::endl;
            } // End of figuring out face data mode


          // read the vertex numbers
          int n_read, v0, v1, v2;
          switch(face_data_mode)
            {
              case 1: n_read = sscanf(line, "%d %d %d", &v0, &v1, &v2); break;
              case 2: n_read = sscanf(line, "%d%*c%*d %d%*c%*d %d%*c%*d", &v0, &v1, &v2); break;
              case 3: n_read = sscanf(line, "%d%*c%*c%*d %d%*c%*c%*d %d%*c%*c%*d ", &v0, &v1, &v2); break;
              case 4: n_read = sscanf(line, "%d%*c%*d%*c%*d %d%*c%*d%*c%*d %d%*c%*d%*c%*d", &v0, &v1, &v2); break;
              default: break;
            }
          if(n_read != 3) printf("mode: %d , line: %s, n_read: %d\n", face_data_mode, line, n_read);
          // printf("%d %d %d \n" , v0, v1, v2);
          assert(n_read == 3);

          // shove the data in the corespoding array
          ftovert.push_back(v0 - 1);
          ftovert.push_back(v1 - 1);
          ftovert.push_back(v2 - 1);
        }
      else if((c == (int)'\n') || (c == (int)' '))
        {
          continue;
        }
      else // if ( c == (int)'#' )
        {
          char line[1000];
          char * success = fgets(line, 1000, fl);
          assert(success == line);
        }
    }

  // close the file
  fclose(fl);

  // Remove the unused vertices
  const int n_uverts = vcoord.size();
  std::vector<bool> is_vert_used(n_uverts, false);

  for(int face = 0; face < int(ftovert.size()); face++)
    {
      is_vert_used[ftovert[face]] = true;
    }

  std::vector<int> uvert_to_vert(n_uverts);
  int n_verts = 0;
  for(int uv = 0; uv < n_uverts; uv++)
    {
      if(is_vert_used[uv])
        {
          uvert_to_vert[uv] = n_verts;
          n_verts++;
        }
      else
        {
          uvert_to_vert[uv] = -1;
        }
    }

  for(int i = 0; i < int(ftovert.size()); i++)
    {
      ftovert[i] = uvert_to_vert[ftovert[i]];
    }

  std::vector<double> vc2(n_verts * 3, -1);
  for(int i = 0; i < n_uverts; i++)
    {
      if(is_vert_used[i])
        {
          vc2[uvert_to_vert[i] * 3 + 0] = vcoord[i * 3 + 0];
          vc2[uvert_to_vert[i] * 3 + 1] = vcoord[i * 3 + 1];
          vc2[uvert_to_vert[i] * 3 + 2] = vcoord[i * 3 + 2];
        }
    }

  // Pass stuff to the grid
  grid().init(vcoord, ftovert, std::vector<int>(), std::vector<int>(), std::vector<int>(), std::vector<int>());
}

//
// Skip the file until you reach a valuable character
//
static bool
skip_nonesense(FILE * fl, const std::string & valuables = "0123456789.-")
{
  bool is_there_more_to_read = false;
  assert(fl);

  for(;;)
    {
      const int c = fgetc(fl);

      if(c == EOF)
        {
          is_there_more_to_read = false;
          break;
        }
      else if(valuables.find((char)c) != std::string::npos)
        {
          int success = fseek(fl, -1, SEEK_CUR);
          is_there_more_to_read = true;
          assert(success == 0);
          break;
        }
    }
  return is_there_more_to_read;
}

//
// Skip the file until you reach a valuable string
//
static std::string
skip_nonesense(FILE * fl, const std::vector<std::string> & tokens)
{

  assert(fl);

  // Put all the valuables in one string and make it unqiue (maybe faster!?)
  std::string valuables;
  {
    for(auto t : tokens)
      {
        valuables += t;
      }
    std::sort(valuables.begin(), valuables.end());
    std::string::iterator last = std::unique(valuables.begin(), valuables.end());
    valuables.erase(last, valuables.end());
  }

  // Skip until you reach a suspect character
  // Then read the whole thing and see if it is what you wanted
  char suspect[1024];
  while(skip_nonesense(fl, valuables))
    {
      fscanf(fl, "%s", suspect);
      for(auto t : tokens)
        {
          if(t == std::string(suspect))
            {
              return t;
            }
        }
    }
  return "";
}


void
Grid_IO::read_off(const std::string & fname)
{

  // open the file
  FILE * fl = open_file(fname, ".off", "r");
  const std::string valuable_characters;

  /*
   * read the file
   */
  int n_verts, n_faces;
  vector<double> vcoord;
  vector<int> ftovert;
  int success;

  // Read the number of verts and faces and edges
  skip_nonesense(fl, valuable_characters);
  success = fscanf(fl, "%d %d %*d", &n_verts, &n_faces);
  assert(success == 2);
  vcoord.reserve(n_verts * 3);
  ftovert.reserve(n_faces * 3);

  // Read vertex coordinates
  for(int v = 0; v < n_verts; v++)
    {
      skip_nonesense(fl, valuable_characters);
      double x, y, z;
      success = fscanf(fl, "%lf %lf %lf", &x, &y, &z);
      assert(success == 3);
      vcoord.push_back(x);
      vcoord.push_back(y);
      vcoord.push_back(z);
    }

  // Read face vertices
  for(int f = 0; f < n_faces; f++)
    {
      skip_nonesense(fl, valuable_characters);
      int nv, v0, v1, v2;
      success = fscanf(fl, "%d %d %d %d", &nv, &v0, &v1, &v2);
      assert(success == 4);
      assert(nv == 3 && "Only triangular faces are supported");
      ftovert.push_back(v0);
      ftovert.push_back(v1);
      ftovert.push_back(v2);
    }

  // close the file
  fclose(fl);

  // init the mesh
  grid().init(vcoord, ftovert, vector<int>(), vector<int>(), vector<int>(), vector<int>());
}

void
Grid_IO::read_gmsh(const std::string & fname)
{
  std::vector<std::string> headers = {"$Nodes", "$Elements"};
  std::string numbers = "0123456789-.";
  int n_nodes, n_elements;
  int node_offset = -1;
  std::vector<double> points;
  std::vector<int> tris;
  std::vector<int> tri_tags;
  std::vector<int> segs;
  std::vector<int> seg_tags;


  FILE * fl = open_file(fname, ".msh", "r");


  //
  // Read the nodes
  //
  io_ensure(skip_nonesense(fl, headers) == "$Nodes", "$Nodes was not found");
  io_ensure(skip_nonesense(fl, numbers), "Unexpected end of file");
  io_ensure(fscanf(fl, "%d", &n_nodes) == 1, "Could not read n_nodes");

  // Read the nodes
  points.reserve(n_nodes);
  for(int i = 0; i < n_nodes; ++i)
    {
      double x, y, z;
      if(node_offset == -1)
        io_ensure(fscanf(fl, "%d", &node_offset) == 1, "Could not read first node index");
      else
        fscanf(fl, "%*d");
      io_ensure(fscanf(fl, "%lf %lf %lf", &x, &y, &z) == 3, "Could not read node coordinates");
      points.push_back(x);
      points.push_back(y);
      points.push_back(z);
    }

  io_ensure(skip_nonesense(fl, std::vector<std::string>({"$EndNodes"})) == "$EndNodes", "$EndNodes was not found");

  //
  // Read the elements
  //
  io_ensure(skip_nonesense(fl, headers) == "$Elements", "$Elements was not found");
  io_ensure(skip_nonesense(fl, numbers), "Unexpected end of file");
  io_ensure(fscanf(fl, "%d", &n_elements) == 1, "Could not read n_nodes");

  // read the elements
  // upper bound reserve
  tris.reserve(n_elements * 3);
  segs.reserve(n_elements * 2);
  tri_tags.reserve(n_elements);
  seg_tags.reserve(n_elements);

  for(int i = 0; i < n_elements; ++i)
    {
      const int element_line = 1;
      const int element_tri = 2;
      const int physical_tag_index = 0;
      const int default_physical_tag = 0;

      int element_type, n_tags, tags[20], physical_tag, v0, v1, v2;

      io_ensure(fscanf(fl, "%*d %d %d", &element_type, &n_tags) == 2, "Could not read elem type and ntags");
      io_ensure(n_tags <= 20, "More tags are not supported");
      for(int j = 0; j < n_tags; ++j)
        {
          io_ensure(fscanf(fl, "%d", tags + j) == 1, "Could not read tag");
        }
      if(n_tags <= physical_tag_index)
        physical_tag = default_physical_tag;
      else
        physical_tag = tags[physical_tag_index];

      if(element_type == element_line)
        {
          io_ensure(fscanf(fl, "%d %d", &v0, &v1) == 2, "Could not read verts 0 and 1");
          seg_tags.push_back(physical_tag);
          segs.push_back(v0 - node_offset);
          segs.push_back(v1 - node_offset);
        }
      else
        {
          io_ensure(element_type == element_tri, "Other element types not supported");
          io_ensure(fscanf(fl, "%d %d %d", &v0, &v1, &v2) == 3, "Could not read verts 0, 1 and 2");
          tri_tags.push_back(physical_tag);
          tris.push_back(v0 - node_offset);
          tris.push_back(v1 - node_offset);
          tris.push_back(v2 - node_offset);
        }
    }

  io_ensure(skip_nonesense(fl, std::vector<std::string>({"$EndElements"})) == "$EndElements",
            "$EndElements was not found");

  // close the file
  fclose(fl);

  // Build the mesh
  grid().init(points, tris, segs, std::vector<int>(), tri_tags, seg_tags);

} // ALl done

void
Grid_IO::read_triangle(const std::string & fname)
{
  const std::string valuables = "0123456789.-";
  int vertex_offset;
  std::vector<int> segs, seg_tags;
  std::vector<int> tris;
  std::vector<double> verts;

  // READ NODES
  {
    int n_vertices, has_tags, n_vertex_attributes;
    FILE * fl = open_file(fname, ".node", "r");

    skip_nonesense(fl, valuables);
    io_ensure(fscanf(fl, "%d %*d %d %d", &n_vertices, &n_vertex_attributes, &has_tags) == 3,
              "Could not read vertex header");

    verts.reserve(3 * n_vertices);
    for(int i = 0; i < n_vertices; ++i)
      {
        int vert_id, tag;
        double x, y, z = 0, attribute;
        skip_nonesense(fl, valuables);
        io_ensure(fscanf(fl, "%d %lf %lf", &vert_id, &x, &y) == 3, "Could not read tri");
        if(i == 0) vertex_offset = vert_id;
        verts.push_back(x);
        verts.push_back(y);
        verts.push_back(z);

        for(int j = 0; j < n_vertex_attributes; ++j)
          {
            io_ensure(fscanf(fl, "%lf", &attribute) == 1, "Could not read attribute");
          }

        if(has_tags)
          {
            io_ensure(fscanf(fl, "%d", &tag) == 1, "Could not read tag");
          }
      }

    fclose(fl);
  }

  // READ VOLUME ELEMENTS
  {
    int n_tris, n_tri_nodes, n_tri_attributes;
    FILE * fl = open_file(fname, ".ele", "r");

    skip_nonesense(fl, valuables);
    io_ensure(fscanf(fl, "%d %d %d", &n_tris, &n_tri_nodes, &n_tri_attributes) == 3, "Could not read vertex header");
    io_ensure(n_tri_nodes == 3, "6 nodes per triangle not supported");

    tris.reserve(3 * n_tris);
    for(int i = 0; i < n_tris; ++i)
      {
        int v0, v1, v2;
        double attribute;
        skip_nonesense(fl, valuables);
        io_ensure(fscanf(fl, "%*d %d %d %d", &v0, &v1, &v2) == 3, "Could not read tri");
        tris.push_back(v0 - vertex_offset);
        tris.push_back(v1 - vertex_offset);
        tris.push_back(v2 - vertex_offset);

        for(int j = 0; j < n_tri_attributes; ++j)
          {
            io_ensure(fscanf(fl, "%lf", &attribute) == 1, "Could not read attribute");
          }
      }

    fclose(fl);
  }

  // READ BOUDNARY ELEMETNS
  {
    int n_vertices, n_segs, has_tags;
    FILE * fl = open_file(fname, ".poly", "r");

    skip_nonesense(fl, valuables);
    io_ensure(fscanf(fl, "%d %*d %*d %*d", &n_vertices) == 1, "Could not read vertex header");
    io_ensure(n_vertices == 0, "Polyfile should not have vertices");

    skip_nonesense(fl, valuables);
    io_ensure(fscanf(fl, "%d %d", &n_segs, &has_tags) == 2, "Could not read segment headers");

    segs.reserve(2 * n_segs);
    if(has_tags) seg_tags.reserve(n_segs);
    for(int i = 0; i < n_segs; ++i)
      {
        int v0, v1, tag;
        skip_nonesense(fl, valuables);
        if(has_tags)
          {
            io_ensure(fscanf(fl, "%*d %d %d %d", &v0, &v1, &tag) == 3, "Could not read segment");
            segs.push_back(v0 - vertex_offset);
            segs.push_back(v1 - vertex_offset);
            seg_tags.push_back(tag);
          }
        else
          {
            io_ensure(fscanf(fl, "%*d %d %d", &v0, &v1) == 2, "Could not read segment");
            segs.push_back(v0 - vertex_offset);
            segs.push_back(v1 - vertex_offset);
          }
      }

    fclose(fl);
  }

  // Create the mesh
  grid().init(verts, tris, segs, std::vector<int>(), std::vector<int>(), seg_tags);
}

void
Grid_IO::read_vtk(const std::string & fl)
{
  assert(0 && "Not supported!!!");
}

void
Grid_IO::write_off(const string & fname) const
{
  Grid_iterators grid_iters(grid());

  // open the file
  FILE * fl = open_file(fname, ".off", "w");

  // write the number of verts and faces and edges
  fprintf(fl, "%d %d %d", grid().n_vertices(), grid().n_real_triangles(), 0);

  // Read vertex coordinates
  for(int v = 0; v < grid_iters.n_vertices(); v++)
    {
      Grid_iterators::Vertex vert = grid_iters.vertex_iterator(v);
      fprintf(fl, "%lf %lf %lf \n", vert.x(), vert.y(), vert.z());
    }
  fprintf(fl, "\n");

  // Read face vertices
  for(int f = 0; f < grid_iters.n_real_triangles(); f++)
    {
      Grid_iterators::Triangle face = grid_iters.triangle_iterator(f);
      const int v1 = face.vertex_from_offset(0).id(); // face.half_edge().origin().id();
      const int v2 = face.vertex_from_offset(1).id(); // face.half_edge().next().origin().id();
      const int v3 = face.vertex_from_offset(2).id(); // face.half_edge().next().next().origin().id();
      fprintf(fl, "%d %d %d %d \n", 3, v1, v2, v3);
    }
  fprintf(fl, "\n");

  // close the file
  fclose(fl);
}

void
Grid_IO::write_vtk(const string & fname)
{
  FILE * fl = open_file(fname, ".vtk", "w");

  //
  this->write_vtk(fl);

  // Write connectivity for debugging
  const int n_rtri = grid().n_real_triangles();
  this->write_vtk_cell_header(fl);
  std::vector<int> v(3 * n_rtri);
  std::vector<int> tag(n_rtri);
  std::vector<int> tagn(3 * n_rtri);
  std::vector<int> n(3 * n_rtri);
  Grid_iterators giter(grid());
  for(int tid = 0; tid < n_rtri; ++tid)
    {
      Grid_iterators::Triangle tri = giter.triangle_iterator(tid);
      v[3 * tid + 0] = tri.vertex_from_offset(0).id();
      v[3 * tid + 1] = tri.vertex_from_offset(1).id();
      v[3 * tid + 2] = tri.vertex_from_offset(2).id();
      tagn[3 * tid + 0] = tri.neighbour_from_offset(0).tag() * (tri.neighbour_from_offset(0).is_dummy() ? -1 : 1);
      tagn[3 * tid + 1] = tri.neighbour_from_offset(1).tag() * (tri.neighbour_from_offset(1).is_dummy() ? -1 : 1);
      tagn[3 * tid + 2] = tri.neighbour_from_offset(2).tag() * (tri.neighbour_from_offset(2).is_dummy() ? -1 : 1);
      n[3 * tid + 0] = tri.neighbour_from_offset(0).raw_id() * (tri.neighbour_from_offset(0).is_dummy() ? -1 : 1);
      n[3 * tid + 1] = tri.neighbour_from_offset(1).raw_id() * (tri.neighbour_from_offset(1).is_dummy() ? -1 : 1);
      n[3 * tid + 2] = tri.neighbour_from_offset(2).raw_id() * (tri.neighbour_from_offset(2).is_dummy() ? -1 : 1);
      tag[tid] = tri.tag();
    }
  this->write_vtk_data(fl, v, "tri_vert_id", true);
  this->write_vtk_data(fl, tagn, "tri_nei_tag", true);
  this->write_vtk_data(fl, n, "tri_nei_id", true);
  this->write_vtk_data(fl, tag, "tri_tag");

  //
  fclose(fl);
}

void
Grid_IO::write_vtk(FILE * fl)
{
  assert(fl);
  Grid_iterators grid_iters(grid());

  /*
   * Write the vtk file.
   */

  // write the header
  fprintf(fl, "# vtk DataFile Version 2.0\n");
  fprintf(fl, "Shayan's output mesh\n");
  fprintf(fl, "ASCII\n");
  fprintf(fl, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fl, "\n");

  // write the vertices
  fprintf(fl, "POINTS %d float\n", grid_iters.n_vertices());
  for(int vidx = 0; vidx < grid_iters.n_vertices(); vidx++)
    {
      Grid_iterators::Vertex vert = grid_iters.vertex_iterator(vidx);
      fprintf(fl, "%e %e %e \n", vert.x(), vert.y(), vert.z());
    }
  fprintf(fl, "\n");

  // write the faces
  fprintf(fl, "CELLS %d %d \n", grid_iters.n_real_triangles(), grid_iters.n_real_triangles() * 4);
  for(int f = 0; f < grid_iters.n_real_triangles(); f++)
    {
      Grid_iterators::Triangle face = grid_iters.triangle_iterator(f);
      int verts[] = {
        face.vertex_from_offset(0).id(),
        face.vertex_from_offset(1).id(),
        face.vertex_from_offset(2).id(),
      };
      fprintf(fl, "3 %d %d %d \n", verts[0], verts[1], verts[2]);
    }
  fprintf(fl, "\n");

  // write the face types
  fprintf(fl, "CELL_TYPES %d \n", grid_iters.n_real_triangles());
  for(int f = 0; f < grid_iters.n_real_triangles(); f++)
    {
      fprintf(fl, "5\n");
    }
  fprintf(fl, "\n");

  _have_written_vert_data_header = _have_written_cell_data_header = false;
}

void
Grid_IO::write_vtk_vert_header(FILE * fl)
{
  if(!_have_written_vert_data_header)
    {
      fprintf(fl, "POINT_DATA %d \n", grid().n_vertices());
      _have_written_vert_data_header = true;
      _writing_mode = MODE_VERTEX;
    }
}

void
Grid_IO::write_vtk_cell_header(FILE * fl)
{
  if(!_have_written_cell_data_header)
    {
      fprintf(fl, "CELL_DATA %d \n", grid().n_real_triangles());
      _have_written_cell_data_header = true;
      _writing_mode = MODE_CELL;
    }
}

void
Grid_IO::write_vtk_data(FILE * fl, std::vector<int> & data, const std::string name, const bool is_vector)
{
  const int bs = (is_vector ? 3 : 1);
  assert(((_writing_mode == (int)MODE_CELL) && ((int)data.size() == grid().n_real_triangles() * bs)) ||
         ((_writing_mode == (int)MODE_VERTEX) && ((int)data.size() == grid().n_vertices() * bs)));

  if(!is_vector)
    {
      fprintf(fl, "SCALARS %s int 1 \n", name.c_str());
      fprintf(fl, "LOOKUP_TABLE default \n");

      for(int i = 0; i < (int)data.size(); i++)
        {
          fprintf(fl, "%d \n", data[i]);
        }
    }
  else
    {
      fprintf(fl, "SCALARS %s int 3 \n", name.c_str());
      fprintf(fl, "LOOKUP_TABLE default \n");

      for(int i = 0; i < (int)(data.size() / 3); i++)
        {
          fprintf(fl, "%d %d %d\n", data[3 * i + 0], data[3 * i + 1], data[3 * i + 2]);
        }
    }
  fprintf(fl, "\n");
}

void
Grid_IO::write_vtk_data(FILE * fl, std::vector<double> & data, const std::string name, const bool is_vector)
{
  const int bs = (is_vector ? 3 : 1);
  assert(((_writing_mode == (int)MODE_CELL) && ((int)data.size() == grid().n_real_triangles() * bs)) ||
         ((_writing_mode == (int)MODE_VERTEX) && ((int)data.size() == grid().n_vertices() * bs)));

  if(!is_vector)
    {
      fprintf(fl, "SCALARS %s float 1 \n", name.c_str());
      fprintf(fl, "LOOKUP_TABLE default \n");

      for(int i = 0; i < (int)data.size(); i++)
        {
          fprintf(fl, "%e \n", data[i]);
        }
    }
  else
    {
      fprintf(fl, "SCALARS %s float 3 \n", name.c_str());
      fprintf(fl, "LOOKUP_TABLE default \n");

      assert(data.size() % 3 == 0);
      for(int i = 0; i < (int)(data.size() / 3); i++)
        {
          fprintf(fl, "%e %e %e\n", data[3 * i + 0], data[3 * i + 1], data[3 * i + 2]);
        }
    }
  fprintf(fl, "\n");
}


void
Grid_IO::print_data_info(FILE * fl) const
{
  Grid_iterators grid_iters(grid());

  {
    typedef std::vector<Grid_connectivity_context::Triangle>::const_iterator Iterator;
    Iterator it = grid().triangles_data().cbegin();
    const Iterator it_end = grid().triangles_data().cend();
    fprintf(fl, "Triangles (%d dummy) \n", grid().n_dummy_triangles());
    for(; it != it_end; ++it)
      {
        fprintf(fl,
                "%d %d %d %d %d %d %d\n",
                it->vert_ids[0],
                it->vert_ids[1],
                it->vert_ids[2],
                it->tri_ids[0],
                it->tri_ids[1],
                it->tri_ids[2],
                it->tag);
      }
  }
  {
    typedef std::vector<Grid_connectivity_context::Vertex>::const_iterator Iterator;
    Iterator it = grid().vertices_data().cbegin();
    Iterator it_end = grid().vertices_data().cend();
    fprintf(fl, "Vertices \n");
    for(; it != it_end; ++it)
      {
        fprintf(fl, "%d %d\n", it->tri_id, it->tag);
      }
  }
}

void
Grid_IO::print_iter_info(FILE * fl) const
{
  assert(fl);
  Grid_iterators grid_iters(grid());

  fprintf(fl, "%8s %8s %8s %8s %8s %8s \n", "#HE", "vbegin", "vend", "face", "twin", "next");

  for(int i = 0; i < grid_iters.n_half_edges(); i++)
    {
      Grid_iterators::Half_edge he = grid_iters.half_edge_iterator(i);
      const int twin = he.twin().raw_half_edge_id();
      const int vbeg = he.origin().id();
      const int vend = he.twin().origin().id();
      const int face = he.triangle().raw_id();
      const int next = he.next().raw_half_edge_id();

      fprintf(fl, "%8d %8d %8d %8d %8d %8d\n", i, vbeg, vend, face, twin, next);
    }
}
