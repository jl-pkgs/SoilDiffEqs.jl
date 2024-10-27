export model_SM_sim, SM_theta2param, of_MSE
export solve_SM_ODE, solve_SM_Bonan

function SM_theta2param(theta)
  n = length(theta) ÷ 2
  κ = theta[1:n]
  cv = theta[n+1:end]
  return κ, cv
end

function model_SM_sim(soil, TS0, theta; method="Bonan", kw...)
  κ, cv = SM_theta2param(theta)
  soil.κ .= κ
  soil.cv .= cv

  if method == "Bonan"
    ysim = solve_SM_Bonan(soil, TS0;)
  elseif method == "ODE"
    ysim = solve_SM_ODE(soil, TS0; kw...)
  end
  ysim
end


function solve_SM_Bonan(soil::Soil{FT}, TS0::AbstractVector{FT}) where {FT<:Real}
  (; n, inds_obs, ibeg) = soil

  ntime = length(TS0)
  R = zeros(ntime, soil.n - ibeg + 1)
  R[1, :] .= soil.θ[ibeg:end]

  for i = 2:ntime
    soil_temperature!(soil, TS0[i]; ibeg)
    R[i, :] .= soil.θ[ibeg:end]
  end

  inds = inds_obs .- ibeg .+ 1
  R[:, inds]
end


"""
    solve_SM_ODE(soil, TS0; solver)

solver = Tsit5()
solver = Rosenbrock23()
solver = Rodas5(autodiff=false)  
"""
function solve_SM_ODE(soil, TS0; solver, reltol=1e-3, abstol=1e-3, verbose=false)
  (; n, inds_obs, ibeg, dt) = soil

  ntime = length(TS0)
  u0 = soil.θ[ibeg:end]

  _SMEquation(dT, T, p, t) = SMEquation_partial(dT, T, p, t; ibeg)
  tspan = (0, dt)
  prob = ODEProblem(_SMEquation, u0, tspan, soil)

  R = zeros(ntime, n - ibeg + 1)
  R[1, :] .= soil.θ[ibeg:end]

  for i = 2:ntime
    soil.TS0 = TS0[i]
    prob.u0 .= soil.θ[ibeg:end]

    sol = solve(prob, solver; reltol, abstol, saveat=dt)
    soil.θ[ibeg:end] .= sol.u[end] # 更新这个时刻的结果
    R[i, :] .= soil.θ[ibeg:end]
  end

  inds = @. inds_obs - ibeg + 1
  R[:, inds]
end
