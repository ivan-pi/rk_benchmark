# Extended notes on callback strategies

Detailed background on each of the six callback idioms benchmarked in this
repository.  See the main [README](README.md) for a quick-reference table and
build/run instructions.

## Strategy 1 — Internal procedure (host association)

The right-hand side is defined as a `contains` subprogram of the calling scope
and captures all needed data through *host association*.  The Fortran 2008
standard explicitly permits passing an internal procedure as an actual
argument, even when the called solver lives in a separate translation unit.

Most compilers implement host-data capture using a **trampoline**: a small,
compiler-generated thunk that holds the host frame pointer and forwards the
call to the real procedure.  Historically these trampolines were written to
the stack at runtime, which required an executable stack — on modern Linux
this triggers a `GNU_STACK` segment marked `RWE` and is rejected by hardened
loaders and some security policies.  Recent GFortran versions (≥ 14) support
`-ftrampoline-impl=heap` to allocate trampolines on the heap instead, and ifx
uses descriptor-based dispatch that avoids the issue entirely.  Whether a
trampoline is actually emitted depends on whether the internal procedure
references host variables; the simple `rhs_internal` in this benchmark does
not capture any locals, so a modern compiler may elide the trampoline
altogether.

The other limitation is structural: the RHS cannot be defined in a separate
module or shared library, which prevents plugin-style extensibility.

## Strategy 2 — RPAR / IPAR arrays (SLATEC style)

Inspired by SLATEC and early NAG/IMSL libraries, the solver passes two
work arrays (`rpar(:)` for reals, `ipar(:)` for integers) through to the
callback unchanged.  The user packs parameters into those arrays.  There is no
dynamic dispatch and no allocation.  The interface is type-safe only in the
sense that the compiler enforces the array types — passing anything other
than real and integer data requires `transfer`, `equivalence`, or similar
type-system subversion, which is cumbersome for structured or heterogeneous
data.

## Strategy 3 — C-style opaque pointer

The callback is declared `bind(C)` and receives a `type(c_ptr)` context
argument (passed by value, C convention).  Inside the callback,
`c_f_pointer` casts the opaque handle back to the concrete derived type.
This is the natural idiom when the solver is a C library and the RHS is
implemented in Fortran, or vice versa.  Dispatch is fully static; the only
overhead is one pointer load via `c_f_pointer`.  This repository optionally
swaps the Fortran implementation for a C one (`rk23_cptr_external.c`,
enabled with `-DENABLE_C_SOLVER=ON`) to demonstrate interoperability in both
directions.

## Strategy 4 — Type-bound procedure / OOP functor

An abstract base type `ode_functor` declares a deferred `eval` method.  The
solver accepts a `class(ode_functor)` argument and calls `fun%eval(...)`.
Problem data is stored as components of the concrete extension type
(`robertson_functor`).  This is the most object-oriented approach: adding a
new problem means extending `ode_functor`, with no changes to the solver
interface.  The cost is dynamic dispatch through the vtable on every call.

## Strategy 5 — Reverse communication interface (RCI)

Instead of calling the user function, the solver *returns* to the caller each
time a function evaluation is needed.  It uses an integer `stage` variable to
encode where it left off; the caller evaluates the RHS, stores the result in
the shared work array, and re-invokes the solver.  There is no callback at
all — the user drives a state machine.  This eliminates all dispatch overhead
and is maximally flexible (the "callback" can be an arbitrary block of code,
including non-Fortran routines or operations that cannot legally appear
inside a callback), but the calling code is verbose and the state-machine
logic must be replicated by every user.  RCI is also the only variant in
this benchmark that naturally supports single-step mode: the caller can
inspect or modify state between any two stages.

## Strategy 6 — Unlimited-polymorphic context `class(*)`

The callback receives a `class(*)` context argument.  Inside the callback a
`select type` construct identifies the dynamic type at runtime and branches
to the appropriate implementation.  This provides a single, fully generic
interface (one abstract interface covers all context types) at the cost of a
runtime type test on every call and the loss of compile-time type safety.

## ODE method

All RK solvers implement the **Bogacki–Shampine RK23** pair with adaptive
step-size control and the FSAL (First Same As Last) property: the fourth
stage evaluation of an accepted step is reused as the first stage of the next
step, saving one function evaluation per accepted step.  Rejected steps do
not consume this saving — `work(:,1)` retains the valid `k1` from the
previous accepted step.

All workspace (`work(neqn, 5)`) is allocated by the caller and passed into
the solver, so no allocation occurs inside the solver loop.  See the main
README for the full list of solver caveats.