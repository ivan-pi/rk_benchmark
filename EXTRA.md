# Extended notes on callback strategies

Detailed background on each of the six callback idioms benchmarked in this
repository.  See the main [README](README.md) for a quick-reference table and
build/run instructions.

## Strategy 1 — Internal procedure (host association)

The right-hand side is defined as a `contains` subprogram of the calling scope
and captures all needed data through *host association*.  However, the compiler
cannot inline the call when the solver lives in a separate translation unit
(compiled library).  In practice, compilers implement the host-data capture
using a *trampoline* — a small, compiler-generated thunk stored on the stack —
which can trigger executable-stack restrictions on some operating systems and
older compiler versions.  The downside is tight coupling: the RHS cannot be
defined in a separate module or shared library.

## Strategy 2 — RPAR / IPAR arrays (SLATEC style)

Inspired by SLATEC and early NAG/IMSL libraries, the solver passes two
work arrays (`rpar(:)` for reals, `ipar(:)` for integers) through to the
callback unchanged.  The user packs parameters into those arrays.  There is no
dynamic dispatch and no allocation.  The interface is type-safe in the sense
that the compiler enforces the array types, but it is constrained to real and
integer data — passing anything else requires `transfer`, `equivalence`, or
other type-system subversion mechanisms, making it cumbersome for complex data.

## Strategy 3 — C-style opaque pointer

The callback is declared `bind(C)` and receives a `type(c_ptr)` context
argument (passed by value, C convention).  Inside the callback,
`c_f_pointer` casts the opaque handle back to the concrete derived type.
This is the natural idiom when the solver is a C library and the RHS is
implemented in Fortran, or vice versa.  It requires no extra dynamic dispatch
(the cast is a simple pointer load) but involves one indirection.

## Strategy 4 — Type-bound procedure / OOP functor

An abstract base type `ode_functor` declares a deferred `eval` method.  The
solver accepts a `class(ode_functor)` argument and calls `fun%eval(...)`.
Problem data is stored as components of the concrete extension type
(`robertson_functor`).  This is the most object-oriented approach: adding a
new problem means extending `ode_functor`, with no changes to the solver
interface.  Dynamic dispatch through the vtable is the cost.

## Strategy 5 — Reverse communication interface (RCI)

Instead of calling the user function, the solver *returns* to the caller each
time a function evaluation is needed.  It uses an integer `stage` variable to
encode where it left off; the caller evaluates the RHS, stores the result in
the shared work array, and re-invokes the solver.  There is **no** callback at
all — the user drives a state machine.  This eliminates all dispatch overhead
and is maximally flexible (the "callback" can be an arbitrary block of code),
but the calling code is verbose and the state machine logic must be replicated
by every user.

## Strategy 6 — Unlimited-polymorphic context `class(*)`

The callback receives a `class(*)` context argument.  Inside the callback a
`select type` construct identifies the dynamic type at runtime and branches to
the appropriate implementation.  This provides a single, fully generic
interface (one abstract interface covers all context types) at the cost of a
runtime branch and the loss of compile-time type safety.

## Repository layout

```
rk_kinds.f90         – kind parameter (dp = c_double) and iso_c_binding imports
rk_types.f90         – abstract interfaces and the ode_functor base type
rk_solvers.f90       – six solver variants (rk23_simple … rk23_class_star)
robertson_models.f90 – Robertson RHS implementations for each strategy
rk_benchmark.f90     – driver: reports mean time per integration over N_runs loops
scripts/             – helper scripts and reference implementations
CMakeLists.txt       – build system
```

## ODE method

All solvers implement the **Bogacki–Shampine RK23** pair with adaptive
step-size control and the FSAL (First Same As Last) property.  The fourth
stage evaluation of an accepted step is reused as the first stage of the next
step, saving one function evaluation per step.

All workspace (`work(neqn, 5)`) is allocated by the caller and passed into the
solver, so no allocation occurs inside the solver loop.
