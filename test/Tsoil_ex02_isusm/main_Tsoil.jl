using Plots, Printf
gr(framestyle=:box)

function plot_soil(i)
  (; z, inds_obs) = soil
  i2 = inds_obs[i]
  title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)
  plot(; title)
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

function plot_obs(i)
  plot(title="layer $i")
  plot!(t, yobs[:, i], label="OBS")
  # plot!(t, TS_sim[:, i], label="SIM")
end


dz = 0.05
nlayer = ceil(Int, 1.35/dz)
Δz = fill(dz, nlayer)
z = [4, 12, 14, 16, 20, 24, 28, 32, 36, 42, 50, 52] * 25.4 / 1000 # inch -> mm -> m
inds_obs = round.(Int, z / dz) # 取这些层的数据

function init_soil(; TS0=20.0, dt=3600.0, soil_type=1, k=3)
  # Δz = fill(0.025, n)
  # Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
  n = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = deepcopy(Tsoil0)

  κ, cv = soil_thermal_properties(Δz, Tsoil, m_liq, m_ice; soil_type)
  Soil{Float64}(; n, ibeg=inds_obs[k], inds_obs=inds_obs[k:end],
    dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, Tsoil)
end

function goal(theta; kw...)
  soil = init_soil(; soil_type=1)
  ysim = model_Tsoil_sim(soil, TS0, theta; kw...)

  obs = yobs[:, 2:end][:]
  sim = ysim[:, 2:end][:]
  # of_MSE(obs, sim)
  gof = GOF(obs, sim)
  -gof.NSE
end
