#include <cstdio>


namespace tests
{
extern bool
point_locator(FILE * fl);
extern bool
connectivity(FILE * fl);
extern bool
numerics(FILE * fl);
extern void
test_poisson(FILE * fl);
extern void
test_transfer(FILE * fl);

const char * const UNIX_COLOR_BLACK = "\033[0;30m";
const char * const UNIX_COLOR_DARKBLUE = "\033[0;34m";
const char * const UNIX_COLOR_DARKGREEN = "\033[0;32m";
const char * const UNIX_COLOR_DARKTEAL = "\033[0;36m";
const char * const UNIX_COLOR_DARKRED = "\033[0;31m";
const char * const UNIX_COLOR_DARKPINK = "\033[0;35m";
const char * const UNIX_COLOR_DARKYELLOW = "\033[0;33m";
const char * const UNIX_COLOR_GRAY = "\033[0;37m";
const char * const UNIX_COLOR_DARKGRAY = "\033[1;30m";
const char * const UNIX_COLOR_BLUE = "\033[1;34m";
const char * const UNIX_COLOR_GREEN = "\033[1;32m";
const char * const UNIX_COLOR_TEAL = "\033[1;36m";
const char * const UNIX_COLOR_RED = "\033[1;31m";
const char * const UNIX_COLOR_PINK = "\033[1;35m";
const char * const UNIX_COLOR_YELLOW = "\033[1;33m";
const char * const UNIX_COLOR_WHITE = "\033[1;37m";
const char * const UNIX_COLOR_RESET = "\033[0m";

static void
print_status(const bool ierr)
{
  if(ierr)
    {
      printf("%s SUCCESS %s \n", UNIX_COLOR_GREEN, UNIX_COLOR_RESET);
    }
  else
    {
      printf("%s FAILED %s \n", UNIX_COLOR_RED, UNIX_COLOR_RESET);
    }
}
}


int
main()
{

  bool ierr;

  printf("=================== Testing connectivity: \n");
  ierr = tests::connectivity(stderr);
  tests::print_status(ierr);
  //
  printf("=================== Testing numerics: \n");
  ierr = tests::numerics(stderr);
  tests::print_status(ierr);
  //
  printf("=================== Testing point locator: \n");
  ierr = tests::point_locator(stderr);
  tests::print_status(ierr);
  //
  // printf("=================== Testing fe_basis lagrange locator: \n");
  // FELagrange::test_lagrange_property();
  // tests::print_status(true);
  //
  printf("=================== Testing Poisson solver: \n");
  tests::test_poisson(stderr);
  tests::print_status(true);

  printf("=================== Testing Grid Transfer: \n");
  tests::test_transfer(stderr);
  tests::print_status(true);



  return 0;
}
