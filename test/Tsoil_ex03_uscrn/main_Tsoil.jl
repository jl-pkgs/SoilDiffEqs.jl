using Plots, Printf
gr(framestyle=:box)

# inner_optimizer = GradientDescent()
# options = Optim.Options(show_trace=true)

z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
Δz = cal_Δz(z)
z, z₊ₕ, dz₊ₕ = soil_depth_init(Δz)
# dz = [2.5, 5, 5, 15, 45, 55]

# Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
# z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

function plot_soil(i; ibeg=2)
  i2 = i + ibeg - 1
  title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)
  plot(; title)
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

function init_soil(; Tsoil0, TS0=20.0, dt=3600.0, soil_type=1, ibeg=2)
  # Δz = fill(0.025, n)
  # Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
  n = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = deepcopy(Tsoil0)

  κ, cv = soil_thermal_properties(Δz, Tsoil, m_liq, m_ice;
    soil_type, method="apparent-heat-capacity")
  Soil{Float64}(; n, dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, Tsoil, ibeg)
end

function goal(theta; yobs, Tsoil0, TS0, ibeg=2, kw...)
  soil = init_soil(; Tsoil0, ibeg, soil_type=7)
  ysim = model_Tsoil_sim(soil, TS0, theta; kw...)

  obs = yobs[:, 2:end]
  sim = ysim[:, 2:end]
  # of_MSE(obs, sim)
  gof = GOF(obs[:], sim[:])
  -gof.NSE
end
