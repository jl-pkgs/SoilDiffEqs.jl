function solve_ode()
  param = ParamVanGenuchten(θs=0.287, θr=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1.0)
  θ0 = 0.267
  ψ0 = van_genuchten_ψ(θ0; param)
  Q0 = -param.Ksat * 0.5  # [cm s-1] 向下为负

  n = 150      # 150 cm?
  Δz = ones(n) # Δz₊ₕ
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹]

  θ0 = fill(0.1, n) |> collect # Example initial soil moisture profile
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation
  p = Soil{Float64}(; n=150, ψ0, θ=θ0, z, z₊ₕ, Δz, Δz₊ₕ, sink, Q0)

  _RichardsEquation(dθ, u, p, t) = RichardsEquation(dθ, u, p, t; method="Q0")
  prob = ODEProblem(_RichardsEquation, θ0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  sol.u[end]
end

function solve_bonan()
  # soil_texture = 1
  param = (soil_texture=1,
    θ_res=0.075, θ_sat=0.287,
    α=0.027, n=3.96, m=1, K_sat=34 / 3600)

  n = 150
  dz = ones(n)
  θ = fill(0.1, n)
  ψ = matric_potential(θ, param; method="van_Genuchten")

  θ0 = 0.267
  ψ0 = matric_potential(θ0, param; method="van_Genuchten")
  Q0 = -param.K_sat * 0.5  # [cm s-1] 向下为负
  
  dt = 5
  ntim = 0.8 * 3600 / dt
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹]
  (; θ, ψ, ψ0, dz, dt, param, sink)

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  # --- Time stepping loop: NTIM iterations with a time step of DT seconds
  for itim = 1:ntim
    hour = itim * (dt / 86400 * 24)
    # @printf("hour = %8.3f\n", hour)
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture_Q0!(θ, ψ, sink, Q0, dz, dt, param)

    # % Sum fluxes for relative mass balance error
    sum_in -= Q0 * dt
    sum_out -= QN * dt
    sum_store += dθ
  end
  SINK = sum(sink) * dt * ntim
  @test (sum_in - sum_out - sum_store - SINK) / sum_in <= 1e4 # 误差小于万分之一
  θ
end

# n = 150
# dz = ones(n)
# z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
@time solution = solve_ode();
@time θ = solve_bonan();

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
