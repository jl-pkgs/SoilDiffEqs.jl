using SoilDifferentialEquations, OrdinaryDiffEq, Test


function data_loader_soil(; TS0=20.0, dt=3600.0, soil_type=1)
  N = 120
  Δz = fill(0.025, N)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = fill(10.0, N)

  κ, cv = soil_properties_thermal(Δz, Tsoil, m_liq, m_ice;
    soil_type, method="apparent-heat-capacity")
  Soil{Float64}(; N, dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, Tsoil)
end


begin
  # 4 hours
  dt = 60 * 6
  soil = data_loader_soil(; dt)
  ntime = round(Int, 3600 * 4 / dt)

  TS_surf = fill(20.0, ntime)
  ysim_bonan = solve_Tsoil_Bonan(soil, TS_surf)

  soil = data_loader_soil(; dt)
  ysim_ode = solve_Tsoil_ODE(soil, TS_surf; solver=Tsit5())

  @test maximum(abs.(ysim_bonan[end, :] - ysim_ode[end, :])) <= 0.05
end

# begin
#   function plot_sm(i)
#     title = "Layer $i"
#     plot(; title)    
#     plot!(ysim_bonan[:, i]; label = "Bonan")
#     plot!(ysim_ode[:, i]; label = "ODE")
#   end
#   layers = [1, 5, 10, 20, 50, 100, 120]
#   plot([plot_sm(i) for i in layers]..., size=(1200, 800))
# end
