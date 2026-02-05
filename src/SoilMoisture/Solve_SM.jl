export ModSim_SM, solve_SM_ODE, solve_SM_Bonan
export SM_param2theta, SM_UpdateParam!, SM_paramBound


function SM_param2theta(soil)
  (; same_layer, method_retention, use_m) = soil.param
  if method_retention == "Campbell"
    (; θ_sat, Ksat, ψ_sat, b) = soil.param
    return same_layer ?
           [θ_sat[1], Ksat[1], ψ_sat[1], b[1]] :
           [θ_sat; Ksat; ψ_sat; b]
  elseif method_retention == "van_Genuchten"
    (; θ_sat, θ_res, Ksat, α, n, m) = soil.param
    theta = same_layer ?
            [θ_sat[1], θ_res[1], Ksat[1], α[1], n[1]] :
            [θ_sat; θ_res; Ksat; α; n]
    m = same_layer ? m[1] : m
    use_m && (theta = [theta; m])
    return theta
  end
end


function SM_UpdateParam!(soil::Soil{T}, theta::AbstractVector{T}) where {T<:Real}
  (; method_retention, same_layer, use_m) = soil.param
  N = soil.N
  if method_retention == "Campbell"
    if same_layer
      soil.param.θ_sat .= theta[1]
      soil.param.Ksat .= theta[2]
      soil.param.ψ_sat .= theta[3]
      soil.param.b .= theta[4]
    else
      soil.param.θ_sat .= theta[1:N]
      soil.param.Ksat .= theta[N+1:2N]
      soil.param.ψ_sat .= theta[2N+1:3N]
      soil.param.b .= theta[3N+1:4N]
    end
  elseif method_retention == "van_Genuchten"
    if same_layer
      soil.param.θ_sat .= theta[1]
      soil.param.θ_res .= theta[2]
      soil.param.Ksat .= theta[3]
      soil.param.α .= theta[4]

      _n = theta[5]
      _m = use_m ? theta[6] : (1 - 1 / _n)
      soil.param.n .= theta[5]
      soil.param.m .= _m
    else
      soil.param.θ_sat .= theta[1:N]
      soil.param.θ_res .= theta[N+1:2N]
      soil.param.Ksat .= theta[2N+1:3N]
      soil.param.α .= theta[3N+1:4N]

      _n = theta[4N+1:5N]
      _m = use_m ? theta[5N+1:6N] : (1 .- 1 ./ _n)
      soil.param.n .= _n
      soil.param.m .= _m
    end
  end
  ## 传递soil.param.param
  Update_SoilParam_Param!(soil.param)
  return nothing
end


function SM_paramBound(soil)
  (; method_retention, same_layer, use_m) = soil.param
  N = soil.N
  if method_retention == "van_Genuchten"
    # θ_sat, θ_res, Ksat, α, n, m
    LOWER = [0.25, 0.03, 0.002, 0.002, 1.05] #, 0.1]
    UPPER = [0.50, 0.20, 60.0, 0.300, 4.00] #, 10.0]
    if use_m
      LOWER = [LOWER; 0.1]
      UPPER = [UPPER; 10.0]
    end
  elseif method_retention == "Campbell"
    # θ_sat, Ksat, ψ_sat, b
    LOWER = [0.25, 0.002, -100, 3.0]
    UPPER = [0.50, 100.0, -5.0, 15.0]
  end

  if same_layer
    return LOWER, UPPER
  else
    return repeat(LOWER; inner=N), repeat(UPPER; inner=N)
  end
end



function ModSim_SM(soil, θ_top; method="Bonan", kw...)
  if method == "Bonan"
    ysim = solve_SM_Bonan(soil, θ_top;)
  elseif method == "ODE"
    ysim = solve_SM_ODE(soil, θ_top; kw...)
  end
  ysim
end


function solve_SM_Bonan(soil::Soil{FT}, θ_top::AbstractVector{FT}) where {FT<:Real}
  (; N, ibeg, inds_obs, sink) = soil
  ntime = length(θ_top)
  R = zeros(ntime, N - ibeg + 1)
  R[1, :] .= soil.θ[ibeg:N]

  for i = 2:ntime
    ψ0 = Init_ψ0(soil, θ_top[i])
    soil_moisture!(soil, sink, ψ0; debug=false)

    @inbounds for j in ibeg:N # copy θ
      j2 = j - ibeg + 1
      R[i, j2] = soil.θ[j]
    end
  end

  inds = @. inds_obs - ibeg + 1
  R[:, inds]
end


"""
    solve_SM_ODE(soil, Tsurf; solver)

solver = Tsit5()
solver = Rosenbrock23()
solver = Rodas5(autodiff=false)  
"""
function solve_SM_ODE(soil, θ_top; solver, reltol=1e-3, abstol=1e-3, verbose=false)
  (; N, inds_obs, ibeg, dt) = soil

  ntime = length(θ_top)
  u0 = soil.θ[ibeg:N]

  _Equation(dθ, θ, p, t) = RichardsEquation_partial(dθ, θ, p, t)
  tspan = (0, dt)
  prob = _ODEProblem(_Equation, u0, tspan, soil)

  R = zeros(ntime, N - ibeg + 1)
  R[1, :] .= soil.θ[ibeg:N]

  for i = 2:ntime
    soil.θ0 = θ_top[i]
    prob.u0 .= soil.θ[ibeg:N]

    sol = _solve(prob, solver; reltol, abstol, saveat=dt)
    soil.θ[ibeg:N] .= sol.u[end] # 更新这个时刻的结果
    R[i, :] .= soil.θ[ibeg:N]
  end

  inds = @. inds_obs - ibeg + 1
  R[:, inds]
end
