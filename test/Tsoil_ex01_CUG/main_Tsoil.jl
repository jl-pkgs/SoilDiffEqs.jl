using Plots, Printf
gr(framestyle=:box)

Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
# _z = [4, 12, 14, 16, 20, 24, 28, 32, 36, 42, 50, 52] * 25.4/1000 # inch to mm

function plot_soil(i; ibeg=1)
  i2 = i + ibeg - 1
  title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)
  plot(; title)
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

function init_soil(; TS0=20.0, dt=3600.0, soil_type=1)
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


# 每次在一个步长上进行求解
function solve_Tsoil_ODE(soil, TS0; reltol=1e-3, abstol=1e-3, verbose=false, ibeg=1)
  ntime = length(TS0)
  R = zeros(ntime, soil.n - ibeg + 1)
  tspan = (0, 3600)
  u0 = soil.Tsoil[ibeg:end]
  _TsoilEquation(dT, T, p, t) = TsoilEquation_partial(dT, T, p, t; ibeg)

  prob = ODEProblem(_TsoilEquation, u0, tspan, soil)
  # prob = ODEProblem(TsoilEquation, u0, tspan, soil)

  # solver = Tsit5()
  # solver = Rosenbrock23()
  solver = Rodas5(autodiff=false)
  R[1, :] .= soil.Tsoil[ibeg:end]

  for i = 2:ntime
    soil.TS0 = TS0[i]
    prob.u0 .= soil.Tsoil[ibeg:end]

    sol = solve(prob, solver; reltol, abstol, saveat=3600)
    soil.Tsoil[ibeg:end] .= sol.u[end] # 更新这个时刻的结果
    R[i, :] .= soil.Tsoil[ibeg:end]
  end
  verbose && (@show soil.timestep)
  R
end


function theta2param(theta)
  n = length(theta) ÷ 2
  κ = theta[1:n]
  cv = theta[n+1:end]
  return κ, cv
end

function model_sim(theta; ibeg::Int=1)
  soil = init_soil(; soil_type=7)
  κ, cv = theta2param(theta)
  soil.κ .= κ
  soil.cv .= cv
  # ysim = solve_Tsoil_ODE(soil, TS0; ibeg)
  ysim = solve_Tsoil_Bonan(soil, TS0; ibeg=1)
  ysim
end

function goal(theta; ibeg::Int=1)
  ysim = model_sim(theta; ibeg)
  obs = yobs[:, 2:end][:]
  sim = ysim[:, 2:end][:]
  # of_MSE(obs, sim)
  gof = GOF(obs, sim)
  -gof.NSE
end

of_MSE(yobs, ysim) = mean((yobs .- ysim) .^ 2)
