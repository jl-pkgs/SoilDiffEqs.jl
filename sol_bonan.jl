using HydroTools, Plots, Test
includet("src/Richards.jl")

function solve_ode()
  param = (; θs=0.287, θr=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1)
  θ0 = 0.267
  ψ0 = van_genuchten_ψ(θ0; param)

  n = 150      # 150 cm?
  Δz = ones(n) # Δz₊ₕ
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹]

  θ0 = fill(0.1, n) |> collect # Example initial soil moisture profile
  tspan = (0.0, 0.8 * 3600)  # Time span for the simulation
  p = Soil{Float64}(; n=150, ψ0, θ=θ0, z, z₊ₕ, Δz, Δz₊ₕ, sink)

  prob = ODEProblem(RichardsEquation, θ0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  sol.u[end]
end
@time solution = solve_ode()

begin
  n = 150
  dz = ones(n)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)

  param = (soil_texture=1,
    θ_res=0.075, θ_sat=0.287,
    α=0.027, n=3.96, m=1, K_sat=34 / 3600)

  θ = fill(0.1, n)
  ψ = matric_potential(θ, param; method="van_Genuchten")

  θ0 = 0.267
  ψ0 = matric_potential(θ0, param; method="van_Genuchten")

  dt = 5
  ntim = 0.8 * 3600 / dt
  sink = ones(n) * 0.3 / 86400 # [cm s⁻¹]
  (; θ, ψ, ψ0, dz, dt, param, sink)

  # % --- Initialize accumulators for water balance check
  sum_in = 0
  sum_out = 0
  sum_store = 0

  # --- Time stepping loop: NTIM iterations with a time step of DT seconds
  @time for itim = 1:ntim
    hour = itim * (dt / 86400 * 24)
    # @printf("hour = %8.3f\n", hour)
    # Calculate soil moisture
    Q0, QN, dθ, err = soil_moisture!(θ, ψ, sink, ψ0, dz, dt, param)

    # % Sum fluxes for relative mass balance error
    sum_in += abs(Q0) * dt
    sum_out += abs(QN) * dt
    sum_store += dθ
  end
  SINK = sum(sink) * dt * ntim
  @test (sum_in - sum_out - sum_store - SINK) / sum_in <= 1e4 # 误差小于万分之一
end

# @test sum_in == 11.810243822643141
# @test sum_out == 0.10508872215771699
# @test sum_store == 11.704825251924781

gr(framestyle=:box)
plot(θ, z; label="Bonan", xlabel="θ", ylabel="Depth (cm)", xlims=(0.08, 0.3))
plot!(solution, z; label="ODE")

maximum(abs.(solution - θ)) <= 1e-3 # 误差小于1/1000
