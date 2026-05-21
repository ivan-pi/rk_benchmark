#!/usr/bin/env python3

"""
Robertson stiff ODE system solved with Diffrax (solver=Bosh3).

This script serves as a Diffrax equivalent to the SciPy solve_ivp (RK23) benchmark.
The problem and tolerances match:
    y0 = [1, 0, 0],  t in [0, 100],  atol = 1e-8,  rtol = 1e-8

Robertson system:
    dy1/dt = -0.04*y1 + 1e4*y2*y3
    dy2/dt =  0.04*y1 - 1e4*y2*y3 - 3e7*y2**2
    dy3/dt =  3e7*y2**2
"""

import jax
import jax.numpy as jnp
from time import perf_counter
from diffrax import diffeqsolve, ODETerm, Bosh3, PIDController, SaveAt

count = {"nfe": 0}

def robertson(t, y, args):
    count["nfe"] += 1
    k1, k2, k3 = 0.04, 1.0e4, 3.0e7
    dy1 = -k1 * y[0] + k2 * y[1] * y[2]
    dy2 =  k1 * y[0] - k2 * y[1] * y[2] - k3 * y[1] ** 2
    dy3 =  k3 * y[1] ** 2
    return jnp.array([dy1, dy2, dy3])


y0 = jnp.array([1.0, 0.0, 0.0])
t0, t1 = 0.0, 100.0
atol = 1.0e-8
rtol = 1.0e-8

term = ODETerm(robertson)
solver = Bosh3()

# PIDController handles the adaptive step sizing to meet atol/rtol requirements
#stepsize_controller = PIDController(rtol=rtol, atol=atol, dcoeff=0)

stepsize_controller = PIDController(
    rtol=rtol, atol=atol,
    pcoeff=0,
    icoeff=1,
    dcoeff=0,
    safety=0.9,
    factormin=0.2,
    factormax=5.0,
)

# Save only the final step. Saving all steps with an explicit solver on a 
# stiff problem results in massive arrays.
saveat = SaveAt(t1=True)

@jax.jit
def run_solve():
    return diffeqsolve(
        term,
        solver,
        t0=t0,
        t1=t1,
        dt0=1.0e-4,
        y0=y0,
        stepsize_controller=stepsize_controller,
        saveat=saveat,
        max_steps=10_000_000  # Explicit solvers need millions of steps for stiff ODEs
    )

print("JIT compiling (warm-up)...")
_ = run_solve()
jax.block_until_ready(_)

# Start benchmark timer
start = perf_counter()
sol = run_solve()
jax.block_until_ready(sol)
elapsed = perf_counter() - start

print("Diffrax diffeqsolve (solver=Bosh3) — Robertson problem [0, 100]")
print(f"  atol = {atol:.0e},  rtol = {rtol:.0e}")
print()
print(f"Elapsed solve time: {elapsed:.6f} s (excluding JIT compilation)")
print()

# Extract the final time and state
final_t = sol.ts[0]
final_y = sol.ys[0]

print(f"Final solution at t = {final_t:.1f}:")
for i, yi in enumerate(final_y):
    print(f"  y[{i+1}] = {yi:.6e}")
print()

# Diffrax tracks step counts in the sol.stats dictionary
accepted_steps = int(sol.stats["num_accepted_steps"])
total_steps = int(sol.stats["num_steps"])
nfev = count["nfe"]

print(f"Accepted steps : {accepted_steps}")
print(f"Total steps    : {total_steps} (including rejected steps)")
print(f"Function eval  : {nfev}")