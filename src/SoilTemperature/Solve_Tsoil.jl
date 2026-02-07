export model_Tsoil_sim, Tsoil_theta2param
export solve_Tsoil_ODE, solve_Tsoil_Bonan
export Tsoil_param2theta, Tsoil_UpdateParam!, Tsoil_paramBound


# theta -> (κ, cv)
function Tsoil_theta2param(theta)
  N = length(theta) ÷ 2
  κ = theta[1:N]
  cv = theta[N+1:end]
  return κ, cv
end

# soil -> theta
function Tsoil_param2theta(soil::Soil{T}) where {T<:Real}
  (; same_layer, κ, cv) = soil.param
  same_layer ? [κ[1], cv[1]] : [κ; cv]
end

function Tsoil_UpdateParam!(soil::Soil{T}, theta::AbstractVector{T}) where {T<:Real}
  N = soil.N
  (; same_layer) = soil.param

  if same_layer
    soil.param.κ .= theta[1]
    soil.param.cv .= theta[2]
  else
    soil.param.κ .= theta[1:N]
    soil.param.cv .= theta[N+1:2N]
  end
end

function Tsoil_paramBound(soil::Soil{T}) where {T<:Real}
  N = soil.N
  (; same_layer) = soil.param

  # κ (热导率): 0.1 - 10 W/m/K
  # cv (热容量): 1e6 - 5e6 J/m³/K
  LOWER = [0.1, 1.0e6]
  UPPER = [10.0, 5.0e6]

  if same_layer
    return LOWER, UPPER
  else
    return repeat(LOWER; inner=N), repeat(UPPER; inner=N)
  end
end


function model_Tsoil_sim(soil, Tsurf, theta; method="Bonan", kw...)
  Tsoil_UpdateParam!(soil, theta)

  if method == "Bonan"
    solve_Tsoil_Bonan(soil, Tsurf)
  else
    solve_Tsoil_ODE(soil, Tsurf; kw...)
  end
end


function solve_Tsoil_Bonan(soil::Soil{FT}, Tsurf::AbstractVector{FT}) where {FT<:Real}
  (; N, inds_obs, ibeg) = soil

  ntime = length(Tsurf)
  R = zeros(ntime, N - ibeg + 1)
  R[1, :] .= soil.Tsoil[ibeg:N]

  for i = 2:ntime
    soil_temperature!(soil, Tsurf[i]; ibeg)
    R[i, :] .= soil.Tsoil[ibeg:N]
  end

  inds = inds_obs .- ibeg .+ 1
  R[:, inds]
end


"""
    solve_Tsoil_ODE(soil, Tsurf; solver, reltol=1e-3, abstol=1e-3)

solver = Tsit5() | Rosenbrock23() | Rodas5(autodiff=false)
"""
function solve_Tsoil_ODE(soil, Tsurf; solver, reltol=1e-3, abstol=1e-3, verbose=false)
  (; N, inds_obs, ibeg, dt) = soil

  ntime = length(Tsurf)
  u0 = soil.Tsoil[ibeg:N]

  _TsoilEquation(dT, T, p, t) = TsoilEquation_partial(dT, T, p, t; ibeg)
  tspan = (0, dt)
  prob = _ODEProblem(_TsoilEquation, u0, tspan, soil)

  R = zeros(ntime, N - ibeg + 1)
  R[1, :] .= soil.Tsoil[ibeg:N]

  for i = 2:ntime
    soil.Tsurf = Tsurf[i]
    prob.u0 .= soil.Tsoil[ibeg:N]

    sol = _solve(prob, solver; reltol, abstol, saveat=dt)
    soil.Tsoil[ibeg:N] .= sol.u[end]
    R[i, :] .= soil.Tsoil[ibeg:N]
  end

  inds = inds_obs .- ibeg .+ 1
  R[:, inds]
end
