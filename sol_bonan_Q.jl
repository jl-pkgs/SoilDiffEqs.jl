using HydroTools, Plots, Test
includet("src/Richards.jl")

function solve_ode()
  param = (; θs=0.287, θr=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1)
  θ0 = 0.267
  ψ0 = van_genuchten_ψ(θ0; param)
  Q0 = -param.Ksat * 0.5  # [cm s-1] 向下为负

  n = 150      # 150 cm?
  Δz = ones(n) # Δz₊ₕ
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹]

  u0 = fill(0.1, n) |> collect # Example initial soil moisture profile
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation
  p = Soil{Float64}(; n=150, ψ0, z, z₊ₕ, Δz, Δz₊ₕ, sink, Q0)

  prob = ODEProblem(RichardsEquation_Q0, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  sol.u[end]
end

n = 150
dz = ones(n)
z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
@time solution = solve_ode()


begin
  # soil_texture = 1
  param = (soil_texture=1,
    Θ_res=0.075, Θ_sat=0.287,
    α=0.027, n=3.96, m=1, K_sat=34 / 3600)

  Θ = fill(0.1, n)
  ψ = matric_potential(Θ, param; method="van_Genuchten")

  Θ0 = 0.267
  ψ0 = matric_potential(Θ0, param; method="van_Genuchten")
  Q0 = -param.K_sat * 0.5  # [cm s-1] 向下为负
  
  dt = 5
  ntim = 0.8 * 3600 / dt
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹]
  (; Θ, ψ, ψ0, dz, dt, param, sink)

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  # --- Time stepping loop: NTIM iterations with a time step of DT seconds
  @time for itim = 1:ntim
    hour = itim * (dt / 86400 * 24)
    # @printf("hour = %8.3f\n", hour)
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture_Q0!(Θ, ψ, sink, Q0, dz, dt, param)

    # % Sum fluxes for relative mass balance error
    sum_in -= Q0 * dt
    sum_out -= QN * dt
    sum_store += dθ
  end
  SINK = sum(sink) * dt * ntim
  @test (sum_in - sum_out - sum_store - SINK) / sum_in <= 1e4 # 误差小于万分之一
end

begin
  gr(framestyle=:box)
  plot(solution, z; label="ODE", xlabel="Θ", ylabel="Depth (cm)", xlims=(0.08, 0.3))
  plot!(Θ, z; label="Bonan")
end

# plot!(solution, z; label="ODE")
maximum(abs.(solution - Θ)) <= 0.003 # 误差小于3/1000
