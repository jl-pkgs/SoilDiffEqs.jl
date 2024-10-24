begin
  n = 120
  Δz = fill(0.025, n)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  dt = 3600
  tspan = (0.0, 3600 * 24)  # Time span for the simulation
end

function init_soil(; TS0=20.0)
  soil_type = 1
  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = fill(2.0 + K0, n)

  κ, cv = soil_thermal_properties(Δz, Tsoil, m_liq, m_ice;
    soil_texture=soil_type, method="apparent-heat-capacity")
  Soil{Float64}(; n=length(z), z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, T=Tsoil)
end

function solve_ode()
  p = init_soil()
  u0 = fill(10.0, p.n)
  p.T .= u0
  prob = ODEProblem(SoilTemperature, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=dt)
  @show p.timestep
  sol.u[end]
end

# solution = "crank-nicolson", "implicit"
function solve_bonan(; TS0=20.0, solution="crank-nicolson")
  soil = init_soil()
  (; Δz, κ, cv) = soil
  Tsoil = fill(10.0, n)

  for t in tspan[1]:dt:tspan[2]
    Tsoil_next, G = soil_temperature(Δz, dt, κ, cv, Tsoil, TS0; solution)
    Tsoil .= Tsoil_next
  end
  Tsoil
end

dt = 3600
tspan = (0.0, dt * 24)  # Time span for the simulation
Tsoil_bonan = solve_bonan()
Tsoil_ode = solve_ode()

maximum(abs.(Tsoil_ode - Tsoil_bonan)) <= 0.05 # 误差小于0.05°

begin
  gr(framestyle=:box)
  plot()
  plot!(Tsoil_bonan, z; label="Bonan")
  plot!(Tsoil_ode, z; label="ODE")
end
