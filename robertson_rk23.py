"""
Robertson stiff ODE system solved with SciPy's solve_ivp (method='RK23').

This script serves as a quick sanity-check reference for the Fortran RK23
benchmark.  The problem and tolerances match those used in rk_benchmark.f90:

    y0 = [1, 0, 0],  t in [0, 100],  atol = 1e-8,  rtol = 1e-8

Robertson system:
    dy1/dt = -0.04*y1 + 1e4*y2*y3
    dy2/dt =  0.04*y1 - 1e4*y2*y3 - 3e7*y2**2
    dy3/dt =  3e7*y2**2
"""

import numpy as np
from scipy.integrate import solve_ivp


def robertson(t, y):
    k1, k2, k3 = 0.04, 1.0e4, 3.0e7
    return [
        -k1 * y[0] + k2 * y[1] * y[2],
         k1 * y[0] - k2 * y[1] * y[2] - k3 * y[1] ** 2,
         k3 * y[1] ** 2,
    ]


y0 = [1.0, 0.0, 0.0]
t_span = (0.0, 100.0)
atol = 1.0e-8
rtol = 1.0e-8

sol = solve_ivp(robertson, t_span, y0, method="RK23", atol=atol, rtol=rtol,
                dense_output=False)

print("SciPy solve_ivp  (method='RK23')  —  Robertson problem  [0, 100]")
print(f"  atol = {atol:.0e},  rtol = {rtol:.0e}")
print()
print(sol)
print()
print(f"Final solution at t = {sol.t[-1]:.1f}:")
for i, yi in enumerate(sol.y[:, -1]):
    print(f"  y[{i+1}] = {yi:.6e}")
print()
print(f"Accepted steps : {sol.t.size - 1}")
print(f"Function evals : {sol.nfev}")
