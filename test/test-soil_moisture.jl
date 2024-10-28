using SoilDifferentialEquations, OrdinaryDiffEq, Test


function data_loader_soil()
  param_water = ParamVanGenuchten(θ_sat=0.287, θ_res=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1.0)
  n = 150
  Δz = fill(0.01, n)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.1, n)
  ψ = van_Genuchten_ψ.(θ; param=param_water)
  θ0 = 0.267
  ψ0 = van_Genuchten_ψ(θ0; param=param_water)

  dt = 5 # [s]
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹], 蒸发速率
  soil = Soil{Float64}(; n, z, z₊ₕ, Δz, Δz₊ₕ, θ, ψ, θ0, ψ0, dt, sink, param_water)
  return soil
end

soil = data_loader_soil()

function solve_ode()
  p = data_loader_soil()
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation
  
  u0 = p.θ
  prob = ODEProblem(RichardsEquation, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  @show p.timestep
  sol.u[end]
end

function solve_bonan()
  soil = data_loader_soil()
  (; dt, ψ0, sink) = soil
  ntim = 0.8 * 3600 / dt

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  # --- Time stepping loop: NTIM iterations with a time step of DT seconds
  for itim = 1:ntim
    hour = itim * (dt / 86400 * 24)
    # @printf("hour = %8.3f\n", hour)
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture!(soil, sink, ψ0)

    # % Sum fluxes for relative mass balance error
    sum_in += abs(Q0) * dt
    sum_out += abs(QN) * dt
    sum_store += dθ
  end
  SINK = sum(sink) * dt * ntim
  @test (sum_in - sum_out - sum_store - SINK) / sum_in <= 1e4 # 误差小于万分之一
  soil.θ
end

@time θ = solve_bonan()
@time solution = solve_ode()

# 40 times slower
# @profview solution = solve_ode();
# @profview for i = 1:20
#   θ = solve_bonan()
# end
@testset "RichardsEquation θ0" begin
  @test maximum(abs.(solution - θ)) <= 0.004 # 误差小于1/1000
end

# @test sum_in == 11.810243822643141
# @test sum_out == 0.10508872215771699
# @test sum_store == 11.704825251924781

# begin
#   using Plots
#   n = 150
#   Δz = fill(0.01, n)
#   z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

#   gr(framestyle=:box)
#   plot(θ, z; label="Bonan", xlabel="θ", ylabel="Depth (cm)", xlims=(0.08, 0.3))
#   plot!(solution, z; label="ODE")
# end
