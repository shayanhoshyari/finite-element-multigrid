cmake_minimum_required(VERSION 3.16)

project(finite_element_multigrid LANGUAGES CXX)

enable_testing()
add_compile_options(-Wall)

set( MG_LIB
  ## Header only
  src/bridson/util.h
  src/bridson/vec.h
  src/bridson/blas_wrappers.h
  src/bridson/sparse_matrix.h # almost a re-write
  src/bridson/pcg_solver.h
  src/bridson/icc.h
  ##
  src/adept/timer.h
  ##
  src/numeric_array.hxx
  src/cell_type.hxx
  src/fake_fprintf.hxx
  src/removable_object.hxx
  ## With .cxx
  src/grid_iterators.hxx
  # src/grid_iterators.cxx INLINED
  src/grid_io.hxx
  src/grid_io.cxx
  src/grid_connectivity_context.hxx
  src/grid_connectivity_context.cxx
  src/point_locator.hxx
  src/point_locator.cxx
  src/quadrature.hxx
  src/quadrature.cxx
  src/poisson_equation.hxx
  src/poisson_equation.cxx
  src/fe.hxx
  src/fe.cxx
  src/mg_tools.hxx
  src/mg_tools.cxx
)
set( MG_LIB_INLINE
  src/grid_iterators.cxx
)

set (CMAKE_CXX_STANDARD 11)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-O0 -g -DDEBUG")
  set(MG_LIB "${MG_LIB}" ${MG_LIB_INLINE})
else()
  add_definitions("-O3 -funroll-loops -march=native -DNDEBUG")
endif()

message("FE_MG_LIB_SOURCES: ${MG_LIB}")

# You probably want to fine tune where this guy finds lapack from.
# On one case it went to my local paraview directories and found
# a random lapack there!!!
# set(LAPACK_DIR "")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
#find_package(CBLAS REQUIRED)
##set(CBLAS_INCLUDE_DIRS "/usr/include/openblas/")
##set(CBLAS_LIBRARIES "/usr/lib/libopenblas.a" -lpthread)
##include_directories(src)
##include_directories("${CBLAS_INCLUDE_DIRS}")

##
## Symlink to meshes
##
function(__add_symlink  __target_name __from __to)
add_custom_target(${__target_name} ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink
                  "${CMAKE_CURRENT_SOURCE_DIR}/${__from}"  "${CMAKE_BINARY_DIR}/${__to}")
endfunction()
__add_symlink(fe_meshes  meshes/dirt meshes)


##
## Executables
##
function(__add_executable  __name )
add_executable( fe_${__name}  exe/${__name}.cxx ${ARGN} ${MG_LIB}  )
#
target_link_libraries(fe_${__name}
${CBLAS_LIBRARIES})
#
 set_target_properties(fe_${__name} PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY         ${CMAKE_BINARY_DIR}/bin
   RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_BINARY_DIR}/bin
   RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin
 )
endfunction()

__add_executable(tests
  exe/test_connectivity.cxx
  exe/test_locator.cxx
  exe/test_numerics.cxx
  exe/test_poisson.cxx
  exe/test_transfer.cxx
  exe/tests.cxx
)

__add_executable(timer_locator
  exe/timer_locator.cxx
)

__add_executable(timer_mg
  exe/timer_mg.cxx
)

__add_executable(benchmark_mg
  exe/benchmark_mg.cxx
)

add_test(NAME finite_element_tests COMMAND fe_tests)
