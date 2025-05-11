import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, perform_step!
  # alg_order, alg_stability_size, explicit_rk_docstring,
  # OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache,
  # alg_cache,
  # OrdinaryDiffEqConstantCache, @fold, trivial_limiter!,
  # constvalue, @unpack, perform_step!, calculate_residuals, @cache,
  # calculate_residuals!, _ode_interpolant, _ode_interpolant!,
  # CompiledFloats, @OnDemandTableauExtract, initialize!,
  # CompositeAlgorithm, _ode_addsteps!, copyat_or_push!,
  # AutoAlgSwitch, get_fsalfirstlast,
  # full_cache, DerivativeOrderNotPossibleError
# using OrdinaryDiffEq, LinearAlgebra

# Define the custom solver type
struct RichardsCNSolver <: OrdinaryDiffEqAlgorithm end

# Core step function for the custom solver
function OrdinaryDiffEq.perform_step!(integrator, cache, alg::RichardsCNSolver)
  # Unpack integrator variables
  u = integrator.u          # Current state (ψ values)
  dt = integrator.dt        # Time step
  p = integrator.p          # Parameters
  t = integrator.t          # Current time
  f = integrator.f          # Derivative function

  # Extract soil parameters and grid info from p
  soil, params, ψ0, sink = p
  N, Δz = soil.N, soil.Δz   # Grid size and spacingspacing

  # Update soil state based on current ψ
  update_soil!(soil, u, params)

  # First step: Predictor (half step)
  A1, b1 = build_tridiagonal_system(soil, dt / 2, ψ0, sink, :predictor)
  u_half = tridiagonal_solve(A1, b1)

  # Update soil state at half step
  update_soil!(soil, u_half, params)

  # Second step: Corrector (full step)
  A2, b2 = build_tridiagonal_system(soil, dt, ψ0, sink, :corrector)
  u_new = tridiagonal_solve(A2, b2)

  # Update the integrator with the new state
  integrator.u = u_new
end

# Placeholder for updating soil properties (to be customized based on your model)
function update_soil!(soil, ψ, params)
  # Update hydraulic conductivity K, water content θ, etc., based on ψ
  # This is model-specific; implement according to your Richards equation setup
end

# Build the tridiagonal system for Crank-Nicolson
function build_tridiagonal_system(soil, dt, ψ0, sink, step_type)
  N = soil.N
  Δz = soil.Δz
  K = soil.K  # Hydraulic conductivity (assumed precomputed in update_soil!)

  # Preallocate tridiagonal matrix components
  a = zeros(N - 1)  # Lower diagonal
  b = ones(N)     # Main diagonal
  c = zeros(N - 1)  # Upper diagonal
  d = copy(ψ0)    # Right-hand side vector

  # Fill tridiagonal system (simplified example)
  for i in 1:N
    if i == 1
      b[i] = 1.0 + (dt / Δz^2) * K[i]
      c[i] = -(dt / Δz^2) * K[i]
      d[i] = ψ0[i] + dt * sink[i]
    elseif i == N
      a[i-1] = -(dt / Δz^2) * K[i]
      b[i] = 1.0 + (dt / Δz^2) * K[i]
      d[i] = ψ0[i] + dt * sink[i]
    else
      a[i-1] = -(dt / Δz^2) * K[i]
      b[i] = 1.0 + 2 * (dt / Δz^2) * K[i]
      c[i] = -(dt / Δz^2) * K[i]
      d[i] = ψ0[i] + dt * sink[i]
    end
  end

  # Return tridiagonal matrix and RHS vector
  A = Tridiagonal(a, b, c)
  return A, d
end

# Thomas algorithm for tridiagonal matrix solving
function tridiagonal_solve(A, b)
end

# Example usage
# Define your ODE function (f!), initial condition (u0), time span (tspan), and parameters (p)
# prob = ODEProblem(f!, u0, tspan, p)
# sol = solve(prob, RichardsCNSolver(), dt=0.01)
