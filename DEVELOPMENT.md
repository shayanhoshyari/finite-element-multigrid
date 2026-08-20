# Development

## 1. Prerequisites

The solver requires a C++11 compiler and CMake 3.16 or newer. Regenerating the meshes requires Jonathan Shewchuk's Triangle 1.6 command-line program. Regenerating the report figures requires Python 3 with Pillow and ReportLab.

## 2. Generate the mesh hierarchy

Place the Triangle executable on `PATH`, or set `TRIANGLE_BIN` explicitly, then run:

```bash
cd meshes
TRIANGLE_BIN=triangle ./generate.sh
cd ..
```

The generated meshes are written to `meshes/dirt/` and intentionally excluded from version control. The eight multigrid meshes contain between 150 and 2,184,494 vertices.

## 3. Build and test

Configure a release build and compile the solver:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build --output-on-failure
```

The benchmark executable is written to `build/bin/fe_benchmark_mg`.

## 4. Reproduce the benchmark

After generating the meshes and building the project, run the full solver matrix:

```bash
python3 report/run_benchmarks.py --build build
```

This runs V-cycle multigrid, W-cycle multigrid, single-grid Gauss–Seidel, GS-preconditioned conjugate gradient, and V-cycle-preconditioned conjugate gradient on mesh levels 2–8. It writes summary data, residual histories, and field samples to `report/generated/`.

The largest cases use more than two million vertices and can take several minutes. In particular, the level-8 single-grid GS run is capped at 2,000 sweeps.

## 5. Regenerate the README figures

Install the plotting dependencies and regenerate the report:

```bash
python3 -m pip install pillow reportlab
python3 report/generate_report.py
```

The generator updates `README.md` and the figures in `report/generated/`.

## 6. Formatting

Format C++ changes with the repository configuration:

```bash
clang-format -i exe/*.cxx src/*.cxx src/*.hxx src/bridson/*.h
```
