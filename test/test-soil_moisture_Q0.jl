using SoilDifferentialEquations, Test
using OrdinaryDiffEqTsit5


function init_soil()
  N = 150
  par = ParamVanGenuchten(θ_sat=0.287, θ_res=0.075, Ksat=34.0, α=0.027, n=3.96, m=1.0)
  param = SoilParam(N, par; use_m=true)

  Δz = fill(0.01, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.1, N)
  ψ = Retention_ψ.(θ; par)
  θ0 = 0.267
  ψ0 = Retention_ψ(θ0, par)
  Q0 = -par.Ksat * 0.5  # [cm h-1] 向下为负

  dt = 5 # [s]
  sink = fill(0.3, N) / 24 # [cm s⁻¹], 蒸发速率
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, θ, ψ,
    Q0, θ0, ψ0, dt, sink, param)
  return soil
end

p = init_soil()
display(p)

function solve_ode()
  p = init_soil()
  u0 = p.θ
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation

  _RichardsEquation(dθ, u, p, t) = RichardsEquation(dθ, u, p, t; method="Q0")
  prob = ODEProblem(_RichardsEquation, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  sol.u[end]
end

function solve_bonan()
  # soil_texture = 1
  soil = init_soil()
  (; dt, Q0, sink) = soil

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  ntim = 0.8 * 3600 / dt
  for itim = 1:ntim
    hour = itim * (dt / 86400 * 24)
    # @printf("hour = %8.3f\N", hour)
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture_Q0!(soil, sink, Q0)

    # % Sum fluxes for relative mass balance error
    sum_in -= Q0 * dt / 3600
    sum_out -= QN * dt / 3600
    sum_store += dθ
  end
  SINK = sum(sink) * dt / 3600 * ntim
  error = sum_in - sum_out - sum_store - SINK
  error_perc = error / sum_in * 100
  @show error, error_perc
  # @test abs(error_perc) <= 1 # 0.2%，千分之二
  soil.θ
end

# N = 150
# dz = ones(N)
# z, z₋ₕ, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
@testset "RichardsEquation Q0" begin
  @time solution = solve_ode()
  @time θ = solve_bonan()
  
  @test maximum(abs.(solution - θ)) <= 2*1e-4 # 误差小于3/1000
end

# begin
#   gr(framestyle=:box)
#   plot(solution, z[1:end]; label="ODE", xlabel="θ", ylabel="Depth (cm)", xlims=(0.08, 0.3))
#   plot!(θ, z[1:end]; label="Bonan")
# end

# plot!(solution, z; label="ODE")
# ODE: 10 times slower
