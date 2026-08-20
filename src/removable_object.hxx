#if !defined(REMOVABLE_OBJECT_HXX_IS_INCLUDED)
#define REMOVABLE_OBJECT_HXX_IS_INCLUDED

#include <cassert>

class Removable_object
{
  public:
  // clang-format off
  Removable_object() : _is_removed(false) {}
  bool is_removed() { return _is_removed; }
  void remove()     { _is_removed = true; }
  void restore()    { assert(_is_removed); _is_removed = false; }
  // clang-format on

private:
  bool _is_removed;
};

#endif /*REMOVABLE_OBJECT_HXX_IS_INCLUDED*/