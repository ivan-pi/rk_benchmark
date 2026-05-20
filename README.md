# rk_benchmark

A Fortran micro-benchmark comparing the runtime cost of different **callback
strategies** inside an adaptive Runge-Kutta ODE solver.

## Objective

Quantify the dispatch overhead of six common Fortran callback idioms.  The
test problem (Robertson chemical kinetics, 3 equations) is intentionally cheap
so that callback overhead is the dominant cost.  Each strategy is timed over
100 complete integrations.

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

## Building and running

**Requirements:** GFortran ≥ 9 (or equivalent Fortran 2008 compiler) and CMake ≥ 3.30.

```bash
git clone https://github.com/ivan-pi/rk_benchmark.git
cd rk_benchmark
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/rk_benchmark
```

To produce a bar-chart PNG of mean time per step:

```bash
gnuplot scripts/plot_mean_time_per_step.gp
# output: scripts/mean_time_per_step.png
```

## Interpreting results

```
    Interface                       Mean(s)   Steps     Rej    NFev     us/step     us/NFev
 1. F77 Ext. (implicit iface)        0.0048      80      59     418     60.3        11.6
 2. Callback with RPAR/IPAR          0.0048      80      59     418     60.4        11.6
 3. Callback C-Style (ctx)           0.0049      80      59     418     61.3        11.7
 4. Functor Method (OOP)             0.0050      80      59     418     62.4        12.0
 5. Reverse Communication            0.0048      80      59     418     60.0        11.5
 6. Class(*) Select Type             0.0051      80      59     418     63.8        12.2
```

* **Mean(s)** — wall time for 100 integrations.
* **Steps / Rej** — accepted and rejected steps (identical across strategies;
  a mismatch would indicate a bug).
* **NFev** — right-hand-side evaluations.
* **us/step** and **us/NFev** — the primary comparison metrics; lower is faster.

For a non-chaotic problem like Robertson kinetics all strategies should produce
identical `Final Y` values (up to floating-point rounding).  Differences
indicate an implementation error.
