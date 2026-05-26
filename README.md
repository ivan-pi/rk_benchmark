# rk_benchmark

A Fortran micro-benchmark comparing the runtime cost of different **callback
strategies** inside an adaptive Runge-Kutta ODE solver.

## Objective

Quantify the dispatch overhead of six common Fortran callback idioms.  The
test problem (Robertson chemical kinetics, 3 equations) is intentionally cheap
so that callback overhead is the dominant cost.

## Tested strategies

See [EXTRA.md](EXTRA.md) for a detailed description of each strategy.

| # | Name | Data passing | Dispatch |
|---|------|-------------|----------|
| 1 | Internal procedure (host association) | Implicit capture (trampoline) | Static |
| 2 | Callback + `rpar`/`ipar` (SLATEC style) | Explicit arrays | Static |
| 3 | Callback + `c_ptr` (C-style opaque pointer) | Opaque pointer | Static (`bind(C)`) |
| 4 | Type-bound procedure (OOP functor) | Object data members | Dynamic (vtable) |
| 5 | Reverse communication interface (RCI) | Caller-managed | None (no callback) |
| 6 | `class(*)` unlimited-polymorphic context | Opaque `class(*)` | Dynamic (`select type`) |

## Methodology caveats

* Results depend heavily on compiler, optimisation flags, and CPU.  Always
  rebuild with `-DCMAKE_BUILD_TYPE=Release` for meaningful numbers.
* The test problem is tiny (3 equations).  Strategies that incur one-time
  setup costs may look worse than they would on larger problems.
* Strategy 1 (internal procedure / trampoline) may trigger executable-stack
  restrictions on some systems.

## Repository layout

```
rk_kinds.f90         – kind parameter (dp = c_double) and iso_c_binding imports
rk_types.f90         – abstract interfaces and the ode_functor base type
rk_solvers.f90       – solver variants plus the direct C binding for the optional external C rk23_cptr solver
rk23_cptr_external.c – optional C implementation of the rk23_cptr solver
robertson_models.f90 – Robertson RHS implementations for each strategy
rk_benchmark.F90     – driver: reports mean time per integration over N_runs loops
rk_forward_euler.f90 – fixed-step Euler runners used by the Euler callback benchmark
euler_benchmark.f90  – driver for the fixed-step forward Euler callback benchmark
scripts/             – helper scripts and reference implementations
CMakeLists.txt       – build system
```

## Building and running

**Requirements:** GFortran ≥ 9 (or equivalent Fortran 2008 compiler), a C compiler, and CMake ≥ 3.30.

```bash
git clone https://github.com/ivan-pi/rk_benchmark.git
cd rk_benchmark
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/rk_benchmark [N_runs]
./build/euler_benchmark [N_runs]
```

To benchmark strategy 3 with the optional external C solver implementation instead
of the Fortran module version, configure with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_C_SOLVER=ON
cmake --build build -j
./build/rk_benchmark [N_runs]
```

To produce a horizontal bar-chart PNG of mean time per step with gnuplot:

```bash
./scripts/plot_mean_time_per_step.gp
# output: scripts/mean_time_per_step.png
```

## Interpreting results

<img src="./results/mean_time_per_step.png" alt="Performance micro-benchmark bar plot" width="600">

```
--------------------------------------------------------------------------------
                         Interface   Mean(s)   Steps     Rej    NFev     us/step
--------------------------------------------------------------------------------
 1. F77 Ext. (implicit iface)         0.0126  137206       8  411643      0.0916
 2. Callback with RPAR/IPAR           0.0127  137206       8  411643      0.0923
 3. Callback C-Style (ctx)            0.0126  137206       8  411643      0.0915
 4. Functor Method (OOP)              0.0127  137206       8  411643      0.0923
 5. Reverse Communication             0.0147  137206       8  411643      0.1071
 6. Class(*) Select Type              0.0125  137206       8  411643      0.0908
--------------------------------------------------------------------------------
```

* **Mean(s)** — mean wall time for a single integration, averaged over N_runs (default: 100) runs.
* **Steps / Rej** — accepted and rejected steps (identical across strategies;
  a mismatch would indicate a bug).
* **NFev** — right-hand-side evaluations.
* **us/step** — the primary comparison metric; lower is faster.

For a non-chaotic problem like Robertson kinetics all strategies should produce
identical `Final Y` values (up to floating-point rounding).  Differences
indicate an implementation error.


## Solver caveats

The six `rk23_*` routines exist to compare callback dispatch overhead, not to serve as production ODE solvers. The following limitations are intentional:

* **Forward integration only.** All variants assume `tend > t` and `h > 0; the termination and endpoint-clipping logic is not direction-aware.
* **I or PI controller only.** Step-size adaptation uses either an I-controller (`fac = 0.9 * err^(-1/3)`) or a PI-controller (`fac = 0.9 * err_n^(-k1/3) * err_{n-1}^(-k2/3)` with defaults `k1=0.7`, `k2=0.4`), both with `[0.2, 5.0]` clamping; no PID controller or stiffness detection.
* **No single-step mode.** Each call integrates from t to tend in one shot; intermediate output requires repeated calls with shorter tend. The reverse-communication variant is the exception by construction.
* **No NaN/Inf checking.** A right-hand side returning non-finite values propagates silently through the error norm and acceptance test.
* **No dense output or event detection.** Bogacki–Shampine admits a cheap cubic Hermite interpolant from the FSAL stages, but it is not exposed; there is no root-finding facility.
* **Explicit method on a mildly stiff problem.** RK23 needs ~137 000 steps for Robertson on `[0, 100]` where a stiff solver would need a few hundred. This is accepted because the goal is callback overhead, not solver quality.