# rk_benchmark

A Fortran micro-benchmark that measures the overhead of different **callback
mechanisms** (strategy-pattern implementations) inside an adaptive
Runge-Kutta ODE solver.

## Table of contents

- [Purpose](#purpose)
- [Motivation](#motivation)
  - [Why the strategy pattern matters](#why-the-strategy-pattern-matters)
  - [Strategy comparison](#strategy-comparison)
    - [Internal procedure (strategy 1)](#internal-procedure-strategy-1)
    - [RPAR / IPAR arrays (strategy 2)](#rpar--ipar-arrays-strategy-2)
    - [C-style opaque pointer (strategy 3)](#c-style-opaque-pointer-strategy-3)
    - [Type-bound procedure / OOP functor (strategy 4)](#type-bound-procedure--oop-functor-strategy-4)
    - [Reverse communication interface (strategy 5)](#reverse-communication-interface-strategy-5)
    - [Unlimited-polymorphic context `class(*)` (strategy 6)](#unlimited-polymorphic-context-class-strategy-6)
- [Repository layout](#repository-layout)
  - [ODE method](#ode-method)
- [Building and running](#building-and-running)

## Purpose

Scientific software frequently needs to decouple a *numerical algorithm* from
the *problem being solved*.  An ODE integrator, for instance, must call a
user-supplied right-hand-side function at every stage of every step, yet it
should know nothing about the physics encoded in that function.  The classical
solution is the **strategy pattern**: the solver receives a callable object (a
*strategy*) that it invokes on demand.

The benchmark isolates the cost of the callback dispatch itself by using an
intentionally cheap test problem — the 3-equation Robertson chemical kinetics
system — so that the integrator overhead represents the largest possible fraction of
total runtime.  Six different strategies are timed over 100 complete
integrations each.

## Motivation

### Why the strategy pattern matters

A numerical library that hard-codes the right-hand side is useless as a
general tool.  The strategy pattern lets the algorithm be written once and
reused across arbitrarily many problems.  In compiled languages this
flexibility comes at a price: every time the solver invokes the user function
it must somehow locate it, and any problem-specific *data* (parameters,
coefficients, lookup tables, …) must travel with the function pointer.

Different languages and standards have arrived at different answers to this
challenge.  Fortran has accumulated several distinct idioms over its fifty-year
history, and each carries different tradeoffs between:

* **Extensibility** — how easy it is to add a new problem without touching
  library code.
* **Data encapsulation** — how naturally problem data travels alongside the
  function.
* **Interoperability** — whether the interface is usable from C or other
  languages.
* **Compiler visibility** — whether the compiler can inline or optimise across
  the callback boundary.
* **Trampoline cost** — whether an extra indirection layer is required.

### Strategy comparison

| # | Name | Data passing | Dispatch | Notes |
|---|------|-------------|----------|-------|
| 1 | Internal procedure (host association) | Implicit capture (trampoline) | Static | Tightest coupling; data lives in the host scope; trampoline may cause executable-stack issues |
| 2 | Callback + `rpar`/`ipar` (SLATEC style) | Explicit arrays | Static | Legacy idiom; constrained to real/integer data; widely portable |
| 3 | Callback + `c_ptr` (C-style opaque pointer) | Opaque pointer | Static (bind(C)) | C-interoperable; requires `c_f_pointer` cast in callback |
| 4 | Type-bound procedure (OOP functor) | Object data members | Dynamic (vtable) | Data and behaviour in one derived type; extensible via inheritance |
| 5 | Reverse communication interface (RCI) | Caller-managed | None (no callback) | Solver yields control; caller evaluates and re-enters; maximally flexible but verbose |
| 6 | `class(*)` unlimited-polymorphic context | Opaque `class(*)` | Dynamic (`select type`) | Single interface for any context type; branch cost inside callback |

#### Internal procedure (strategy 1)

The right-hand side is defined as a `contains` subprogram of the calling scope
and captures all needed data through *host association*.  However, the compiler
cannot inline the call when the solver lives in a separate translation unit
(compiled library).  In practice, compilers implement the host-data capture
using a *trampoline* — a small, compiler-generated thunk stored on the stack —
which can trigger executable-stack restrictions on some operating systems and
older compiler versions.  The downside is tight coupling: the RHS cannot be
defined in a separate module or shared library.

#### RPAR / IPAR arrays (strategy 2)

Inspired by SLATEC and early NAG/IMSL libraries, the solver passes two
work arrays (`rpar(:)` for reals, `ipar(:)` for integers) through to the
callback unchanged.  The user packs parameters into those arrays.  There is no
dynamic dispatch and no allocation.  The interface is type-safe in the sense
that the compiler enforces the array types, but it is constrained to real and
integer data — passing anything else requires `transfer`, `equivalence`, or
other type-system subversion mechanisms, making it cumbersome for complex data.

#### C-style opaque pointer (strategy 3)

The callback is declared `bind(C)` and receives a `type(c_ptr)` context
argument (passed by value, C convention).  Inside the callback,
`c_f_pointer` casts the opaque handle back to the concrete derived type.
This is the natural idiom when the solver is a C library and the RHS is
implemented in Fortran, or vice versa.  It requires no extra dynamic dispatch
(the cast is a simple pointer load) but involves one indirection.

#### Type-bound procedure / OOP functor (strategy 4)

An abstract base type `ode_functor` declares a deferred `eval` method.  The
solver accepts a `class(ode_functor)` argument and calls `fun%eval(...)`.
Problem data is stored as components of the concrete extension type
(`robertson_functor`).  This is the most object-oriented approach: adding a
new problem means extending `ode_functor`, with no changes to the solver
interface.  Dynamic dispatch through the vtable is the cost.

#### Reverse communication interface (strategy 5)

Instead of calling the user function, the solver *returns* to the caller each
time a function evaluation is needed.  It uses an integer `stage` variable to
encode where it left off; the caller evaluates the RHS, stores the result in
the shared work array, and re-invokes the solver.  There is **no** callback at
all — the user drives a state machine.  This eliminates all dispatch overhead
and is maximally flexible (the "callback" can be an arbitrary block of code),
but the calling code is verbose and the state machine logic must be replicated
by every user.

#### Unlimited-polymorphic context `class(*)` (strategy 6)

The callback receives a `class(*)` context argument.  Inside the callback a
`select type` construct identifies the dynamic type at runtime and branches to
the appropriate implementation.  This provides a single, fully generic
interface (one abstract interface covers all context types) at the cost of a
runtime branch and the loss of compile-time type safety.

## Repository layout

```
rk_kinds.f90        – kind parameter (dp = c_double) and iso_c_binding imports
rk_types.f90        – abstract interfaces and the ode_functor base type
rk_solvers.f90      – six solver variants (rk23_simple … rk23_class_star)
robertson_models.f90– Robertson RHS implementations for each strategy
rk_benchmark.f90    – driver program: times each variant over N_runs loops
CMakeLists.txt      – build system
```

### ODE method

All solvers implement the **Bogacki–Shampine RK23** pair with adaptive
step-size control and the FSAL (First Same As Last) property.  The fourth
stage evaluation of an accepted step is reused as the first stage of the next
step, saving one function evaluation per step.

All workspace (`work(neqn, 5)`) is allocated by the caller and passed into the
solver, so no allocation occurs inside the solver loop.

## Building and running

**Requirements:** a Fortran 2008 (or later) compiler such as GFortran ≥ 9,
and CMake ≥ 3.30.

```bash
git clone https://github.com/ivan-pi/rk_benchmark.git
cd rk_benchmark
cmake -B build -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build
./build/rk_benchmark
```

## Plot mean time per step

From the repository root, run:

```bash
gnuplot scripts/plot_mean_time_per_step.gp
```

This runs `./build/rk_benchmark`, pipes the benchmark table through `awk`, and
writes a bar chart to `scripts/mean_time_per_step.png`.

Example output (timings will vary by machine):

```
RK23 Final Refactored Benchmark (Clean FSAL Property)
Integrations per test: 100
--------------------------------------------------------------------------------
1. Internal Proc (Host-Data):      0.4821 s | Final Y:      1.0000      0.0000      0.0000
2. Callback with RPAR/IPAR:        0.4835 s | Final Y:      1.0000      0.0000      0.0000
3. Callback C-Style (ctx):         0.4902 s | Final Y:      1.0000      0.0000      0.0000
4. Functor Method (OOP):           0.5013 s | Final Y:      1.0000      0.0000      0.0000
5. Reverse Communication:          0.4798 s | Final Y:      1.0000      0.0000      0.0000
6. Class(*) Select Type:           0.5142 s | Final Y:      1.0000      0.0000      0.0000
```
