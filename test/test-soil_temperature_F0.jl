using SoilDifferentialEquations, OrdinaryDiffEq, Test


function init_soil(; TS0=20.0, dt=3600.0, soil_type=1)
  N = 120
  Δz = fill(0.025, N)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = fill(10.0, N)

  F0 = -10.0 # [W m-2]
  κ, cv = soil_properties_thermal(Δz, Tsoil, m_liq, m_ice;
    soil_type, method="apparent-heat-capacity")
  Soil{Float64}(; N, dt, z, z₊ₕ, Δz, Δz₊ₕ, TS0, Tsoil, F0, 
    param = SoilParam(; N, κ, cv))
end

function solve_ode(reltol=1e-5, abstol=1e-5)
  p = init_soil()
  (; dt, Tsoil) = p
  u0 = Tsoil
  _TsoilEquation(dT, T, p, t) = TsoilEquation(dT, T, p, t; method="F0")

  prob = ODEProblem(_TsoilEquation, u0, tspan, p)
  sol = solve(prob, Tsit5(); reltol, abstol, saveat=dt)
  @show p.timestep
  sol.u[end]
end

# solution = "crank-nicolson", "implicit"
function solve_bonan(; solution="crank-nicolson")
  soil = init_soil()
  (; dt, Tsoil) = soil

  ts = tspan[1]:dt:tspan[2]
  ntime = length(ts)
  df0 = 0.0
  f0 = -soil.F0
  length(f0) == 1 && (f0 = fill(f0, ntime))

  for k in 1:ntime
    soil_temperature_F0!(soil, df0, f0[k]) # Tsoil_next, G = 
  end
  Tsoil
end

soil = init_soil()
tspan = (0.0, 3600 * 24 * 10)  # Time span for the simulation

@testset "Tsoil" begin
  @time Tsoil_bonan = solve_bonan()
  @time Tsoil_ode = solve_ode()
  @test maximum(abs.(Tsoil_ode - Tsoil_bonan)) <= 0.007 # 0.0131
end

# begin
#   gr(framestyle=:box)
#   plot(xlabel="Temperature [°C]", ylabel="Depth [cm]", yflip=true)
#   plot!(Tsoil_bonan, -z*100; label="Bonan")
#   plot!(Tsoil_ode, -z*100; label="ODE")
# end
