#ifndef FAKE_FPRINTF_IS_INCLUDED
#define FAKE_FPRINTF_IS_INCLUDED

#include <cstdarg>
#include <cstdio>
#include <iostream>

typedef int (*Fprintf_function)(FILE *, const char *, ...);

inline int
fake_fprintf(FILE *, const char *, ...)
{
  return 0;
}


struct File_stream
: std::ostream
, std::streambuf
{
  File_stream(FILE * fl)
  : std::ostream(this)
  , _fl(fl)
  {
    if(!_fl)
      {
        printf("FILE must be open");
        exit(-1);
      }
  }

  int overflow(int c) { return putc(c, _fl); }

private:
  FILE * _fl;
};

#endif /*FAKE_FPRINTF_IS_INCLUDED*/