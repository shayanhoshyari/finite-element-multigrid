#ifndef CELL_TYPE_IS_INCLUDED
#define CELL_TYPE_IS_INCLUDED

enum CellType
{
  CELL_TYPE_TET = 0,
  CELL_TYPE_HEX,
  CELL_TYPE_WEDGE,
  CELL_TYPE_PYRAMID,
  CELL_TYPE_TRI,
  CELL_TYPE_QUAD,
  CELL_TYPE_EDGE
};

inline int
cell_dim(const CellType type)
{
  switch(type)
    {
      case CELL_TYPE_WEDGE:
      case CELL_TYPE_PYRAMID:
      case CELL_TYPE_TET:
      case CELL_TYPE_HEX: return 3;
      case CELL_TYPE_TRI:
      case CELL_TYPE_QUAD: return 2;
      case CELL_TYPE_EDGE: return 1;
    }

  throw;
  return -1;
}

inline int
cell_n_verts(const CellType type)
{
  switch(type)
    {
      case CELL_TYPE_WEDGE: return 6;
      case CELL_TYPE_PYRAMID: return 5;
      case CELL_TYPE_TET: return 4;
      case CELL_TYPE_HEX: return 8; ;
      case CELL_TYPE_TRI: return 3;
      case CELL_TYPE_QUAD: return 4;
      case CELL_TYPE_EDGE: return 2;
    }

  throw;
  return -1;
}

inline const char *
cell_type_name(const CellType type)
{
  static const char * __CELL_TYPE_WEDGE = "WEDGE";
  static const char * __CELL_TYPE_PYRAMID = "PYRAMID";
  static const char * __CELL_TYPE_TET = "TET";
  static const char * __CELL_TYPE_HEX = "HEX";
  static const char * __CELL_TYPE_TRI = "TRI";
  static const char * __CELL_TYPE_QUAD = "QUAD";
  static const char * __CELL_TYPE_EDGE = "EDGE";

  switch(type)
    {
      case CELL_TYPE_WEDGE: return __CELL_TYPE_WEDGE;
      case CELL_TYPE_PYRAMID: return __CELL_TYPE_PYRAMID;
      case CELL_TYPE_TET: return __CELL_TYPE_TET;
      case CELL_TYPE_HEX: return __CELL_TYPE_HEX;
      case CELL_TYPE_TRI: return __CELL_TYPE_TRI;
      case CELL_TYPE_QUAD: return __CELL_TYPE_QUAD;
      case CELL_TYPE_EDGE: return __CELL_TYPE_EDGE;
    }

  throw;
  return "";
}


#define N_MAX_CELL_NODES 8

#endif /*CELL_TYPE_IS_INCLUDED*/