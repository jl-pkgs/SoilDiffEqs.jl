## kongdd, 2025-02-14
# 1. dt只有时间步长比较小，才能取得较高的精度
using SoilDifferentialEquations, Test
using OrdinaryDiffEqTsit5
# using Plots

function init_soil()
  N = 150
  # par = VanGenuchten(θ_sat=0.287, θ_res=0.075, Ksat=34.0, α=0.027, n=3.96, m=1.0)
  par = VanGenuchten(θ_sat=0.287, θ_res=0.075, Ksat=14.0, α=0.027, n=3.96, m=1.0)
  # par = VanGenuchten(θ_sat=0.387, θ_res=0.075, Ksat=34.0, α=0.075, n=1.8)
  param = SoilParam(N, par; use_m=true)

  Δz = fill(0.01, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.1, N)
  ψ = Retention_ψ.(θ; par) # cm
  θ0 = 0.267
  ψ0 = Retention_ψ(θ0; par)
  
  dt = 30 # seconds
  sink = fill(0.3, N) / 24 # [cm h⁻¹], 蒸发速率0.3cm/d
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, θ, ψ, θ0, ψ0, dt, sink, param)
  return soil
end

function solve_ode()
  p = init_soil()
  tspan = (0.0, 0.8*3600)  # Time span for the simulation
  u0 = p.θ
  prob = ODEProblem(RichardsEquation, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  @show p.timestep
  return sol.u[end]
end

function solve_bonan()
  soil = init_soil()
  (; dt, ψ0, sink) = soil
  ntim = 0.8 * 3600 / dt

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  # --- Time stepping loop: NTIM iterations with a time step of DT seconds
  for itim = 1:ntim
    # hour = itim * dt / 3600 # in hour
    # @printf("hour = %8.3f\N", hour)
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture!(soil, sink, ψ0)
    # println("Q0 = $Q0, QN = $QN, dθ = $dθ, err = $err")
    
    # % Sum fluxes for relative mass balance error
    sum_in += abs(Q0) * dt / 3600
    sum_out += abs(QN) * dt / 3600
    sum_store += dθ
  end
  
  SINK = sum(sink) * dt / 3600 * ntim
  error = sum_in - sum_out - sum_store - SINK
  error_perc = error / sum_in * 100
  @show error, error_perc
  @test abs(error_perc) <= 1 # 0.2%，千分之二
  return soil.θ
end

@testset "soil_moisture!" begin
  @time θ = solve_bonan()
  @time solution = solve_ode()
  @test maximum(abs.(solution - θ)) * 100 <= 0.2 # 误差小于0.2%
end

# begin
#   using Plots
#   N = 150
#   Δz = fill(0.01, N)
#   z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

#   gr(framestyle=:box)
#   plot(θ, z[1:end]; label="Bonan", xlabel="θ", ylabel="Depth (cm)", xlims=(0.08, 0.3))
#   plot!(solution, z[1:end]; label="ODE")
# end

# 40 times slower
# @profview solution = solve_ode();
# @profview for i = 1:20
#   θ = solve_bonan()
# end
# @testset "RichardsEquation θ0" begin
#   
# end

# @test sum_in == 11.810243822643141
# @test sum_out == 0.10508872215771699
# @test sum_store == 11.704825251924781
