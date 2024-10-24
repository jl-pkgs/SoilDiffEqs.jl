gr(framestyle=:box)

function plot_soil(i)
  plot(title="layer $i")
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

function init_soil(; TS0=20.0, dt=3600.0, soil_type=1)
  # Δz = fill(0.025, n)
  Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
  n = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = fill(10.0, n)

  κ, cv = soil_thermal_properties(Δz, Tsoil, m_liq, m_ice;
    soil_texture=soil_type, method="apparent-heat-capacity")
  # κ, cv：两个比较重要的参数
  Soil{Float64}(; n, dt, z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, TS0, Tsoil)
end


# 每次在一个步长上进行求解
function solve_Tsoil_ODE(soil, TS0; reltol=1e-5, abstol=1e-5, verbose=false)
  ntime = length(TS0)
  R = zeros(ntime, soil.n)
  # soil = init_soil()
  tspan = (0, 3600)
  u0 = soil.Tsoil
  prob = ODEProblem(TsoilEquation, u0, tspan, soil)

  for i = 1:ntime
    soil.TS0 = TS0[i]
    prob.u0 .= soil.Tsoil
    # @assert prob.p.TS0 == TS0[i] # 确认
    sol = solve(prob, Tsit5(); reltol, abstol, saveat=3600)
    soil.Tsoil .= sol.u[end] # 更新这个时刻的结果
    R[i, :] .= soil.Tsoil
  end
  verbose && (@show soil.timestep)
  # soil.Tsoil
  R
end


function theta2param(theta)
  κ = fill(theta[1], 9)
  cv = fill(theta[2], 9)
  if length(theta) == 2
    κ = fill(theta[1], 9)
    cv = fill(theta[2], 9)  
  end
  return κ, cv
end

function model_sim(theta)
  soil = init_soil(; soil_type=7)
  κ, cv = theta2param(theta)
  soil.κ .= κ
  soil.cv .= cv
  ysim = solve_Tsoil_ODE(soil, TS0)
  ysim
end

function goal(theta)
  ysim = model_sim(theta)
  obs = yobs[:, 2:end][:]
  sim = ysim[:, 2:end][:]
  # of_MSE(obs, sim)
  gof = GOF(obs, sim)
  -gof.NSE
end

of_MSE(yobs, ysim) = mean((yobs .- ysim) .^ 2)
