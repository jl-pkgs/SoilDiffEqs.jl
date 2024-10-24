function init_soil(; TS0=20.0, dt=3600.0, soil_type=1)
  n = 120
  Δz = fill(0.025, n)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = fill(10.0, n)

  κ, cv = soil_thermal_properties(Δz, Tsoil, m_liq, m_ice;
    soil_texture=soil_type, method="apparent-heat-capacity")
  Soil{Float64}(; n, dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, Tsoil)
end


function solve_Tsoil_ODE(TS0; reltol=1e-5, abstol=1e-5)
  soil = init_soil()
  # 每次在一个步长上进行求解
  tspan = (0, 3600)
  u0 = soil.Tsoil
  prob = ODEProblem(TsoilEquation, u0, tspan, soil)

  ntime = length(TS0)
  for i = 1:ntime
    soil.TS0 = TS0[i]
    prob.u0 .= soil.Tsoil
    # @assert prob.p.TS0 == TS0[i] # 确认
    sol = solve(prob, Tsit5(); reltol, abstol, saveat=3600)
    soil.Tsoil .= sol.u[end]
  end
  @show soil.timestep
  soil.Tsoil
end

# solution = "crank-nicolson", "implicit"
function solve_bonan(TS0; solution="crank-nicolson")
  soil = init_soil()
  (; dt, Δz, κ, cv, Tsoil) = soil

  ntime = length(TS0)
  for k in 1:ntime
    Tsoil_next, G = soil_temperature(Δz, dt, κ, cv, Tsoil, TS0[k]; solution)
    Tsoil .= Tsoil_next
  end
  Tsoil
end

# TS0 = fill(23.0, nrow(d))
TS0 = A[:, 1]
@time Tsoil_bonan = solve_bonan(TS0);
@time Tsoil_ode = solve_Tsoil_ODE(TS0);
@test maximum(abs.(Tsoil_ode - Tsoil_bonan)) <= 0.35

# @profview Tsoil_ode = solve_ode()
begin
  n = 120
  Δz = fill(0.025, n)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  gr(framestyle=:box)
  plot(xlabel="Temperature [°C]", ylabel="Depth [cm]", yflip=true)
  plot!(Tsoil_bonan, -z * 100; label="Bonan")
  plot!(Tsoil_ode, -z * 100; label="ODE")
end
