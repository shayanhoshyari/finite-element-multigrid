# Source API Guide

## Contents

- [1. Scope](#1-scope)
- [2. Component map](#2-component-map)
- [3. Mesh storage](#3-mesh-storage)
- [4. Mesh I/O](#4-mesh-io)
- [5. Point location and mesh transfer](#5-point-location-and-mesh-transfer)
- [6. Finite elements and quadrature](#6-finite-elements-and-quadrature)
- [7. Poisson problem and assembly](#7-poisson-problem-and-assembly)
- [8. Numeric and sparse types](#8-numeric-and-sparse-types)
- [9. PCG and preconditioners](#9-pcg-and-preconditioners)
- [10. Multigrid](#10-multigrid)
- [11. Complete minimal workflow](#11-complete-minimal-workflow)
- [12. Known API limitations](#12-known-api-limitations)

## 1. Scope

The `src/` directory is a small C++11 finite-element and geometric-multigrid
library for scalar Poisson problems on two-dimensional triangular meshes. Its
main responsibilities are:

1. storing triangle-mesh connectivity;
2. navigating that connectivity through vertices, triangles, and half-edges;
3. reading and writing meshes;
4. evaluating finite-element basis functions and quadrature;
5. assembling a sparse Poisson system;
6. locating one mesh in another and transferring nodal fields; and
7. applying multigrid as either a solver or a PCG preconditioner.

The API is research-oriented. Most contracts are enforced with `assert`, so a
Debug build is recommended while developing. Some public-looking methods are
unfinished; these are called out below.

## 2. Component map

| Component | Primary types | Responsibility |
|---|---|---|
| Mesh storage | `Grid_connectivity_context` | Own vertices, triangles, adjacency, boundary tags, and boundary dummy triangles |
| Mesh navigation | `Grid_iterators` | Lightweight vertex, triangle, half-edge, and vertex-neighbour views |
| Mesh I/O | `Grid_IO` | Read OBJ, OFF, Gmsh, VTK, or Triangle meshes; write OFF or VTK |
| Geometry search | `Point_locator` | Locate points, project to boundaries, and locate the vertices of one mesh in another |
| Quadrature and basis data | `QGenerator`, `FE_segment`, `FE_triangle` | Reference-cell quadrature and physical-cell linear basis data |
| PDE definition and assembly | `Poisson_problem_data`, `Poisson_assembler` | Describe a scalar Poisson problem and assemble its sparse system |
| Numeric storage | `Numeric_array`, `bridson::FixedSparseMatrix` | Fixed-size-after-initialization vectors and CSR-like sparse matrices |
| Iterative methods | `bridson::PCGSolver`, `bridson::Preconditioner` | PCG and reusable preconditioners |
| Multigrid | `mg::Parameters`, `mg::Context` | Mesh hierarchy, transfer, smoothing, cycles, and solve history |

The usual dependency flow is:

```text
mesh files -> Grid_connectivity_context -> Poisson_assembler -> A and f
                          |
                          +-> Point_locator -> inter-grid transfer data

A, f + transfer data + smoothers -> mg::Context -> solution
```

## 3. Mesh storage

### 3.1. `Grid_connectivity_context`

Declared in `src/grid_connectivity_context.hxx`.

`Grid_connectivity_context` owns a triangular mesh. A real triangle stores
three vertex IDs and three neighbouring-triangle IDs. Each boundary edge is
represented internally by an additional dummy triangle; this lets the
half-edge navigation use the same topology at interior and boundary edges.

The principal construction API is:

```cpp
Grid_connectivity_context grid;

Grid_connectivity_context grid(
    const std::vector<double>& vertices,
    const std::vector<int>& triangles,
    const std::vector<int>& boundary_segments = {},
    const std::vector<int>& vertex_tags = {},
    const std::vector<int>& triangle_tags = {},
    const std::vector<int>& boundary_segment_tags = {});
```

`vertices` is flattened in blocks of three coordinates, `triangles` in blocks
of three vertex IDs, and `boundary_segments` in blocks of two vertex IDs.

Useful queries include:

```cpp
int n_vertices() const;
int n_real_triangles() const;
int n_dummy_triangles() const;
int n_all_triangles() const;

const std::vector<Vertex>& vertices_data() const;
const std::vector<Triangle>& triangles_data() const;
const std::vector<int>& boundary_tags() const;
```

The type supports deep copying through its copy constructor and `set_equal`,
and efficient content exchange through `swap`.

`reorder_vertices()` is not implemented and must not be called.

### 3.2. `Grid_iterators`

Declared in `src/grid_iterators.hxx`.

`Grid_iterators` is a non-owning navigation view over a
`Grid_connectivity_context`:

```cpp
Grid_iterators iterators(grid);
auto vertex = iterators.vertex_iterator(vertex_id);
auto triangle = iterators.triangle_iterator(triangle_id);
auto edge = iterators.boundary_half_edge_iterator(boundary_edge_id);
```

Its nested view types are inexpensive handles containing a pointer to the mesh
and one or two integer IDs:

- `Vertex` exposes `id`, `tag`, coordinates, an incident triangle or half-edge,
  and boundary queries.
- `Triangle` exposes vertices, neighbours, its first half-edge, its tag, and
  whether it is a boundary dummy triangle.
- `Half_edge` exposes `origin`, `next`, `prev`, `twin`, and owning triangle.
- `Vertex_umbrella` walks the neighbours incident on a vertex.

A typical neighbour traversal is:

```cpp
Grid_iterators::Vertex_umbrella umbrella =
    iterators.vertex_umbrella_iterator(vertex_id);
do {
  Grid_iterators::Vertex neighbour = umbrella.half_edge().origin();
  // use neighbour
} while (umbrella.advance());
```

These are navigation handles, not STL iterators: they do not define
`operator++`, ranges, or ownership. Combinatorial mesh editing is not
implemented.

## 4. Mesh I/O

`Grid_IO`, declared in `src/grid_io.hxx`, wraps a mutable mesh reference:

```cpp
Grid_connectivity_context grid;
Grid_IO io(grid);
io.read_triangle("meshes/squaremg.4");
io.write_vtk("mesh.vtk");
```

The supported readers are selected explicitly with `read_obj`, `read_off`,
`read_gmsh`, `read_vtk`, and `read_triangle`, or by filename through
`read_auto`. Writers are provided for OFF and legacy VTK.

VTK field output is stateful: write the mesh, then the appropriate vertex or
cell header, and then one or more arrays with `write_vtk_data`.

```cpp
FILE* output = std::fopen("solution.vtk", "w");
io.write_vtk(output);
io.write_vtk_vert_header(output);
io.write_vtk_data(output, solution, "solution");
std::fclose(output);
```

I/O failures and malformed data are primarily handled through assertions or
process-level errors rather than exceptions.

## 5. Point location and mesh transfer

### 5.1. `Point_locator`

Declared in `src/point_locator.hxx`.

A locator is a non-owning view over the mesh to be searched:

```cpp
Point_locator locator(search_grid);
```

The API has three levels:

- static or mesh-bound triangle and edge searches;
- point queries using directional, breadth-first, or brute-force search; and
- whole-mesh searches used to build multigrid transfer information.

Queries return one of `LOCATOR_SUCCESS`, `LOCATOR_FAILURE`, or
`LOCATOR_FATAL_ERROR`. `Tri_query_out` contains the owner triangle and
barycentric coordinates. `Edge_query_out` contains the boundary edge,
one-dimensional coordinate, and squared distance.

For mesh transfer, the important call is:

```cpp
Point_locator::Grid_query_out transfer;
int status = Point_locator(coarse_grid).bfs_mesh_search(fine_grid, transfer);
```

This locates every fine-grid vertex in the coarse grid. `transfer.owner_id[i]`
identifies its coarse owner triangle, while the two entries beginning at
`transfer.xi[2*i]` contain its coordinates. A boundary owner may be a dummy
triangle, in which case the coordinate convention is edge-specific.

The breadth-first whole-mesh search is the intended fast path. The brute-force
version is most useful as a reference implementation in tests.

### 5.2. Transfer operations

`mg::Context` exposes transfer between adjacent hierarchy levels:

```cpp
context.prolongate_dofs(coarse_level, fine_level, coarse_values, fine_values);
context.restrict_dofs(fine_level, coarse_level, fine_values, coarse_values);
```

Levels are numbered from finest (`0`) to coarsest (`n_levels - 1`). The two
levels must be adjacent. Prolongation is barycentric interpolation, and
restriction is its transpose.

The static `interpolate_impl` and `interpolate_transpose_impl` functions can be
used without constructing a complete multigrid hierarchy, but their argument
direction must match the `Grid_query_out` convention described above.

## 6. Finite elements and quadrature

### 6.1. `QGenerator`

`QGenerator`, declared in `src/quadrature.hxx`, constructs quadrature points
and weights for a `CellType` and polynomial order:

```cpp
QGenerator quadrature;
quadrature.init(CELL_TYPE_TRI, 3);

for (unsigned i = 0; i < quadrature.n_points(); ++i) {
  const bridson::Vec3d& point = quadrature.qp(i);
  double weight = quadrature.w(i);
}
```

The generator contains rules for edges, triangles, quadrilaterals,
tetrahedra, prisms, pyramids, and hexahedra. The rest of this project primarily
uses edges and triangles.

### 6.2. `FE_segment` and `FE_triangle`

Declared in `src/fe.hxx`.

These classes represent linear Lagrange basis functions evaluated at a fixed
quadrature order. Construct once, then call `reinit` for each physical cell:

```cpp
FE_triangle fe(3);
fe.reinit(p0, p1, p2);

const auto& physical_points = fe.qp();
const auto& weighted_jacobians = fe.JxW();
const auto& basis_values = fe.phi();
const auto& basis_gradients = fe.dphi();
```

Array indexing is `phi[basis_id][quadrature_point_id]`, with the same layout
for `dphi`. `FE_segment` provides basis values but does not populate basis
gradients.

## 7. Poisson problem and assembly

### 7.1. `Poisson_problem_data`

Declared in `src/poisson_equation.hxx`.

Users describe a scalar problem by deriving from `Poisson_problem_data`:

```cpp
class Problem : public Poisson_problem_data {
public:
  bridson::Vec4d stiffness() const override;
  double right_hand_side(const bridson::Vec2d& x) const override;
  void boundary_condition_terms(
      const bridson::Vec2d& x,
      int boundary_tag,
      double& q,
      double& t) const override;
};
```

The assembled equation is

```text
-div(K grad(u)) = f
```

with Robin boundary data

```text
q(x) u + grad(u) . n = t(x).
```

`stiffness()` returns the flattened 2-by-2 tensor `K`. Exact-solution support
is optional: override `has_exact_solution()` and `exact_solution()` when it is
available.

The problem object must outlive any `Poisson_assembler` referring to it.

### 7.2. `Poisson_assembler`

The assembler determines the nodal sparsity structure and fills the matrix and
right-hand side:

```cpp
Problem problem;
Poisson_assembler assembler(problem);

bridson::FixedSparseMatrix A;
assembler.get_nonzero_structure(grid, A);

Numeric_array rhs(assembler.get_n_dof(grid));
assembler.assemble_system(grid, A, rhs);
```

There is one degree of freedom per vertex. Calling
`get_nonzero_structure` before `assemble_system` is required. If an exact
solution is provided, `find_exact_dofs` evaluates it at every mesh vertex.

## 8. Numeric and sparse types

### 8.1. `Numeric_array`

`Numeric_array` is a public subclass of `std::vector<double>` intended to have
a fixed size after initialization. Construct it with a size, or default
construct it and call one of the `init` methods exactly once.

It provides common solver operations:

```cpp
x.set_zero();
x.set_equal(y);
x.add(y);
x.subtract(y);
x.add_scaled(alpha, y);

double l1 = x.norm_l1();
double l2 = x.norm_l2();
double linf = x.norm_linf();
```

The fixed-size contract is assertion-based and can still be bypassed through
inherited `std::vector` mutators. `Numeric_array::dot(other)` currently ignores
`other`; use `bridson::BLAS::dot(x, y)` for a cross-vector dot product.

### 8.2. `bridson::FixedSparseMatrix`

This is a square CSR-like matrix with public storage:

- `xadj[i]` and `xadj[i+1]` delimit row `i`;
- `adj` contains column IDs; and
- `values` contains the corresponding entries.

Call `set_adjacency` before setting or adding entries. It sorts column IDs
within each row, which is required by the Gauss-Seidel preconditioner.

```cpp
bridson::FixedSparseMatrix A;
A.set_adjacency(row_offsets, column_ids);
A.set_zero();
A.add_values(row_ids, column_ids, local_matrix);
```

Free functions `multiply` and `multiply_scale_and_add` provide sparse
matrix-vector products.

## 9. PCG and preconditioners

`bridson::Preconditioner` defines the interface:

```cpp
struct Preconditioner {
  virtual Preconditioner* clone() const = 0;
  virtual void form(const FixedSparseMatrix& matrix) = 0;
  virtual void apply(const std::vector<double>& rhs,
                     std::vector<double>& result) = 0;
};
```

Available implementations are:

- `Dummy_preconditioner`;
- `Jacobi_preconditioner`;
- `ICC_preconditioner`; and
- `Gauss_seidel_preconditioner` (symmetric Gauss-Seidel).

PCG usage is:

```cpp
bridson::PCGSolver solver(
    new bridson::Gauss_seidel_preconditioner,
    true); // solver takes ownership

solver.set_solver_parameters(relative_tolerance, max_iterations);

double residual;
int iterations;
bool converged = solver.solve(A, rhs, solution, residual, iterations);
```

The tolerance is relative to the initial infinity norm of the residual, and
`residual` is also an infinity norm. `solver.history` retains the initial and
per-iteration residuals.

The pointer constructor only owns the preconditioner when its second argument
is `true`. Prefer a lifetime that makes this ownership explicit.

## 10. Multigrid

### 10.1. Hierarchy conventions

`mg::Context`, declared in `src/mg_tools.hxx`, owns the complete multigrid
state. Levels are ordered as follows:

```text
level 0                         finest grid and target system
level 1
...
level n_levels - 1             coarsest grid
```

Its public vectors expose the grids, system matrices, smoothers, work arrays,
and transfer metadata. This is useful for experiments, although callers must
preserve their size and level ordering.

### 10.2. Required setup order

The intended lifecycle is:

```cpp
mg::Parameters parameters;
parameters.cycle_type = mg::Parameters::CYCLE_TYPE_V;
parameters.relaxation_factor = 1.0;
parameters.n_per_level_smoothing = {3, 3, 500};

mg::Context context;
context.set_parameters(parameters);
context.set_grids(grids);                  // swaps from finest to coarsest
context.set_up_discrete_equations(assembler);
context.set_up_transfer_info();
context.set_smoothers(bridson::Gauss_seidel_preconditioner());
context.verify();
```

`set_grids` swaps the contents out of the input vector. The input must be
ordered finest to coarsest, and every successive grid must have fewer
vertices. `set_smoothers(const Preconditioner&)` clones and forms one smoother
per level. Its pointer-vector overload takes ownership of every pointer and
sets the input entries to null.

`relaxation_factor` is currently stored but not applied by the implementation.

### 10.3. Cycles and standalone solve

One V- or W-cycle is applied with:

```cpp
double previous_residual = context.cycle(level, mg::ZERO);
```

Most clients should instead use the solver interface:

```cpp
context.set_solver_parameters(1e-10, 40);

double residual;
int iterations;
bool converged = context.solve(rhs, solution, residual, iterations);
```

The multigrid stopping tolerance is an absolute Euclidean residual norm.
`context.history` contains one residual value per completed cycle.

At each non-coarse level, a cycle performs pre-smoothing, residual
restriction, one recursive call for a V-cycle or two for a W-cycle,
prolongation and correction, and post-smoothing. The coarsest level applies
the configured number of smoother steps.

### 10.4. Multigrid as a PCG preconditioner

`mg::Context` implements `bridson::Preconditioner`, so a configured hierarchy
can be passed to PCG without transferring ownership:

```cpp
bridson::PCGSolver solver(&context, false);
solver.set_solver_parameters(1e-10, 100);
solver.solve(context.AA[0], rhs, solution, residual, iterations);
```

Each `apply` starts from a zero correction and performs one multigrid cycle.
`mg::Context::clone()` is not implemented, so the overload of
`set_smoothers` that clones a preconditioner cannot be used with another
`mg::Context`.

## 11. Complete minimal workflow

The benchmark in `exe/benchmark_mg.cxx` is the canonical integration example.
In condensed form, a standalone solve looks like this:

```cpp
Problem problem;
Poisson_assembler assembler(problem);

std::vector<Grid_connectivity_context> grids(3);
Grid_IO(grids[0]).read_triangle("meshes/squaremg.3");
Grid_IO(grids[1]).read_triangle("meshes/squaremg.2");
Grid_IO(grids[2]).read_triangle("meshes/squaremg.1");

mg::Parameters parameters;
parameters.cycle_type = mg::Parameters::CYCLE_TYPE_V;
parameters.relaxation_factor = 1.0;
parameters.n_per_level_smoothing = {3, 3, 500};

mg::Context context;
context.set_parameters(parameters);
context.set_grids(grids);
context.set_up_discrete_equations(assembler);
context.set_up_transfer_info();
context.set_smoothers(bridson::Gauss_seidel_preconditioner());

Numeric_array rhs(context.ff[0][mg::ZERO].size());
Numeric_array solution(rhs.size());
rhs.set_equal(context.ff[0][mg::ZERO]);

context.set_solver_parameters(1e-10, 40);
double residual;
int iterations;
bool converged = context.solve(rhs, solution, residual, iterations);
```

## 12. Known API limitations

- Assertions carry many precondition and error checks and disappear in
  optimized builds.
- `Grid_connectivity_context::reorder_vertices()` is not implemented.
- Mesh combinatorial editing is not supported.
- `mg::Context::clone()` is not implemented.
- `mg::Parameters::relaxation_factor` is not used.
- `Numeric_array::dot(other)` does not currently use `other`.
- The point-locator boundary coordinate convention differs from its interior
  barycentric convention.
- Mesh I/O is C-style and does not provide structured error objects.
- The library does not currently expose a stable installed target, namespace
  all public types, or promise ABI compatibility.

These constraints are manageable for the repository's current role as a
finite-element multigrid study, but callers should not treat the source API as
a hardened general-purpose library interface.
