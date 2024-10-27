export model_SM_sim, SM_theta2param, of_MSE
export solve_SM_ODE, solve_SM_Bonan


function model_SM_sim(soil, θ_surf, theta; method="Bonan", kw...)
  soil.K .= theta

  if method == "Bonan"
    ysim = solve_SM_Bonan(soil, θ_surf;)
  elseif method == "ODE"
    ysim = solve_SM_ODE(soil, θ_surf; kw...)
  end
  ysim
end


function solve_SM_Bonan(soil::Soil{FT}, θ_surf::AbstractVector{FT}) where {FT<:Real}
  (; n, inds_obs, ibeg) = soil

  ntime = length(θ_surf)
  R = zeros(ntime, soil.n - ibeg + 1)
  R[1, :] .= soil.θ[ibeg:end]

  for i = 2:ntime
    soil_temperature!(soil, θ_surf[i]; ibeg)
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
function solve_SM_ODE(soil, θ_surf; solver, reltol=1e-3, abstol=1e-3, verbose=false)
  (; n, inds_obs, ibeg, dt) = soil

  ntime = length(θ_surf)
  u0 = soil.θ[ibeg:end]

  _Equation(dθ, θ, p, t) = RichardsEquation(dθ, θ, p, t; ibeg)
  tspan = (0, dt)
  prob = ODEProblem(_Equation, u0, tspan, soil)

  R = zeros(ntime, n - ibeg + 1)
  R[1, :] .= soil.θ[ibeg:end]

  for i = 2:ntime
    soil.θ0 = θ_surf[i]
    prob.u0 .= soil.θ[ibeg:end]

    sol = solve(prob, solver; reltol, abstol, saveat=dt)
    soil.θ[ibeg:end] .= sol.u[end] # 更新这个时刻的结果
    R[i, :] .= soil.θ[ibeg:end]
  end

  inds = @. inds_obs - ibeg + 1
  R[:, inds]
end
