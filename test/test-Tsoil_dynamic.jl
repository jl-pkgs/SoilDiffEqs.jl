function init_soil(; Tsurf=20.0, dt=3600.0, soil_type=1)
  N = 120
  Δz = fill(0.025, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = fill(10.0, N)

  κ, cv = soil_properties_thermal(Δz, Tsoil, m_liq, m_ice;
    soil_texture=soil_type, method="apparent-heat-capacity")
  Soil{Float64}(; N, dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, Tsurf, Tsoil)
end


function solve_Tsoil_ODE(Tsurf; reltol=1e-5, abstol=1e-5)
  soil = init_soil()
  # 每次在一个步长上进行求解
  tspan = (0, 3600)
  u0 = soil.Tsoil
  prob = ODEProblem(TsoilEquation, u0, tspan, soil)

  ntime = length(Tsurf)
  for i = 1:ntime
    soil.Tsurf = Tsurf[i]
    prob.u0 .= soil.Tsoil
    # @assert prob.p.Tsurf == Tsurf[i] # 确认
    sol = solve(prob, Tsit5(); reltol, abstol, saveat=3600)
    soil.Tsoil .= sol.u[end] # 更新这个时刻的结果
  end
  @show soil.timestep
  soil.Tsoil
end

# solution = "crank-nicolson", "implicit"
function solve_bonan(Tsurf; solution="crank-nicolson")
  soil = init_soil()
  (; dt, Δz, κ, cv, Tsoil) = soil

  ntime = length(Tsurf)
  for k in 1:ntime
    Tsoil_next, G = soil_temperature(Δz, dt, κ, cv, Tsoil, Tsurf[k]; solution)
    Tsoil .= Tsoil_next
  end
  Tsoil
end

# Tsurf = fill(23.0, nrow(d))
Tsurf = A[:, 1]
@time Tsoil_bonan = solve_bonan(Tsurf);
@time Tsoil_ode = solve_Tsoil_ODE(Tsurf);
@test maximum(abs.(Tsoil_ode - Tsoil_bonan)) <= 0.35

# @profview Tsoil_ode = solve_ode()
begin
  N = 120
  Δz = fill(0.025, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  gr(framestyle=:box)
  plot(xlabel="Temperature [°C]", ylabel="Depth [cm]", yflip=true)
  plot!(Tsoil_bonan, -z * 100; label="Bonan")
  plot!(Tsoil_ode, -z * 100; label="ODE")
end
