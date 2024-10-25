function init_soil(Δz; TS0=20.0, dt=3600.0, soil_type=1)
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
  # κ, cv：两个比较重要的参数
  Soil{Float64}(; n, dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, Tsoil)
end

function theta2param(theta)
  n = length(theta) ÷ 2
  κ = fill(theta[1], n)
  cv = fill(theta[2], n)
  return κ, cv
end

param2theta(soil) = vcat(soil.κ, soil.cv)

of_MSE(yobs, ysim) = mean((yobs .- ysim) .^ 2)
