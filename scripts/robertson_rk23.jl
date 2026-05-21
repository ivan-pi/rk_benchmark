using DifferentialEquations
using Printf

# Define the Robertson system
function robertson!(du, u, p, t)
    k1, k2, k3 = 0.04, 1.0e4, 3.0e7

    du[1] = -k1 * u[1] + k2 * u[2] * u[3]
    du[2] =  k1 * u[1] - k2 * u[2] * u[3] - k3 * u[2]^2
    du[3] =  k3 * u[2]^2
end

# Parameters
y0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
atol = 1.0e-8
rtol = 1.0e-8

# Setup and solve
prob = ODEProblem(robertson!, y0, tspan)

# Warmup run to compile
solve(prob, BS3(), abstol=atol, reltol=rtol);

# BS3 is the standard Bogacki-Shampine implementation in DiffEq.jl
# The default controller is PI (Proportional-Integral) for this solver.
start_time = time()
sol = solve(prob, BS3(), abstol=atol, reltol=rtol)
elapsed = time() - start_time

# Output results
@printf("Julia DifferentialEquations.jl (method='BS3') — Robertson problem [0, 100]\n")
@printf("  atol = %.0e,  rtol = %.0e\n\n", atol, rtol)
@printf("Elapsed solve time: %.6f s\n\n", elapsed)
@printf("Final solution at t = %.1f:\n", sol.t[end])
for i in 1:3
    @printf("  y[%d] = %.6e\n", i, sol.u[end][i])
end
@printf("\nAccepted steps : %d\n", length(sol.t) - 1)
@printf("Function evals : %d\n", sol.destats.nf)
