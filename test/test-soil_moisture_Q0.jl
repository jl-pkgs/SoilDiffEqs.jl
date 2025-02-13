using SoilDifferentialEquations, Test
using OrdinaryDiffEq


function data_loader_soil()
  N = 150
  par = ParamVanGenuchten(θ_sat=0.287, θ_res=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1.0)
  param = SoilParam(N, par; use_m=true)

  Δz = fill(0.01, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.1, N)
  ψ = van_Genuchten_ψ.(θ; par)
  θ0 = 0.267
  ψ0 = van_Genuchten_ψ(θ0, par)
  Q0 = -par.Ksat * 0.5  # [cm s-1] 向下为负

  dt = 5 # [s]
  sink = ones(N) * 0.3 / 86400 # [cm s⁻¹], 蒸发速率
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, θ, ψ,
    Q0, θ0, ψ0, dt, sink, param)
  return soil
end

p = data_loader_soil()
display(p)

function solve_ode()
  p = data_loader_soil()
  u0 = p.θ
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation

  _RichardsEquation(dθ, u, p, t) = RichardsEquation(dθ, u, p, t; method="Q0")
  prob = ODEProblem(_RichardsEquation, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  sol.u[end]
end

function solve_bonan()
  # soil_texture = 1
  soil = data_loader_soil()
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
    sum_in -= Q0 * dt
    sum_out -= QN * dt
    sum_store += dθ
  end
  SINK = sum(sink) * dt * ntim
  @test (sum_in - sum_out - sum_store - SINK) / sum_in <= 1e4 # 误差小于万分之一
  soil.θ
end

# N = 150
# dz = ones(N)
# z, z₋ₕ, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
@time solution = solve_ode()
@time θ = solve_bonan()

@testset "RichardsEquation Q0" begin
  @test maximum(abs.(solution - θ)) <= 0.003 # 误差小于3/1000
end

# begin
#   gr(framestyle=:box)
#   plot(solution, z; label="ODE", xlabel="θ", ylabel="Depth (cm)", xlims=(0.08, 0.3))
#   plot!(θ, z; label="Bonan")
# end
# plot!(solution, z; label="ODE")
# ODE: 10 times slower
