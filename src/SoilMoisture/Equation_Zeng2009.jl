"""
- Q0: Infiltration rate, [cm h-1], should be negative or zero
"""
function cal_Q_Zeng2009!(soil::Soil{T}, θ::AbstractVector{T}) where {T<:Real}
  (; N, jwt, Q, ψ, ψE, K₊ₕ) = soil
  (; θ_sat, θ_res, param) = soil.param
  zwt = soil.zwt * 100 # [m] -> [cm]
  z = soil.z_cm
  Δz = soil.Δz_cm

  cal_θEψE!(soil) # update θE, ψE
  cal_K!(soil, θ)
  cal_ψ!(soil, θ)

  z[N+1] = 0.5 * (zwt + z[N])
  Δz[N+1] = zwt >= z[N] ? Δz[N] : abs(zwt - z[N])

  i = N
  if jwt == N # GW under soil profile
    _θ = 0.5 * (θ[i] + θ_sat[i])
    se = clamp((_θ - θ_res[i]) / (θ_sat[i] - θ_res[i]), 0.01, 1.0)
    # se = 0.5 * (θ_sat[i] + θ[i]) / θ_sat[i]
    # se = clamp(se, 0.01, 1.0)
    ψ[i+1] = Retention_ψ_Se(se, param[i]) # N+1层的ψ，用的是第N层
  end
  jwt < N && (ψ[i+1] = 0.0)

  ## kernal
  for i = 1:N
    dz = (z[i+1] - z[i])
    dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
    Q[i] = -K₊ₕ[i] * dψ / dz
  end
  jwt < N && (Q[N] = 0.0)
  return Q
end


function RichardsEquation_Zeng2009(dθ::AbstractVector{T}, θ::AbstractVector{T}, soil::Soil{T}, t; method="ψ0") where {T<:Real}
  # 这里采用的是Q0
  soil.timestep += 1
  (; ibeg, N, Q, Q0, sink) = soil # Δz, z, 
  Δz = soil.Δz_cm
  cal_Q_Zeng2009!(soil, θ)

  dθ[ibeg] = ((-Q0 + Q[ibeg]) - sink[ibeg]) / Δz[ibeg] / 3600.0
  @inbounds for i in ibeg+1:N
    dθ[i] = ((-Q[i-1] + Q[i]) - sink[i]) / Δz[i] / 3600.0 # [m3 m-3] / h-1
  end
end


"""
    solve_SM_ODE(soil, Tsurf; solver)

solver = Tsit5()
solver = Rosenbrock23()
solver = Rodas5(autodiff=false)  
"""
function solve_SM_Zeng2009(soil; solver, reltol=1e-3, abstol=1e-3, verbose=false,
  ET::AbstractVector)

  ntime = length(ET)
  (; N, inds_obs, ibeg, dt) = soil

  ## 加入蒸发分配模块, 蒸发总量
  u0 = soil.θ[ibeg:N]

  # _Equation(dθ, θ, p, t) = RichardsEquation_Zeng2009(dθ, θ, p, t)
  tspan = (0, dt)
  prob = _ODEProblem(RichardsEquation_Zeng2009, u0, tspan, soil)

  SM = zeros(ntime, N - ibeg + 1)
  SINK = zeros(ntime, N)
  Q = zeros(ntime, N)
  SM[1, :] .= soil.θ[ibeg:N]

  nchunk = 50
  chunksize = ceil(Int, ntime / nchunk)
  nchunk = ceil(Int, ntime / chunksize)
  p = Progress(nchunk)

  # 加入土壤含水量限制因子
  for i = 2:ntime
    (verbose && mod(i, chunksize) == 0) && next!(p)
    # soil.θ0 = θ_surf[i]
    prob.u0 .= soil.θ[ibeg:N]

    k = 2
    β = soil.θ[k] / soil.param.θ_sat[k]
    soil.sink[k] = ET[i] / 10 * β # mm/h to cm/h

    soil.θ_prev[1:N] .= soil.θ[1:N]
    sol = _solve(prob, solver; reltol, abstol, saveat=dt)
    soil.θ[ibeg:N] .= sol.u[end] # 更新这个时刻的结果

    SINK[i, :] .= soil.sink[1:N]
    Q[i, :] = cal_Q_Zeng2009!(soil, soil.θ)
    SM[i, :] .= soil.θ[ibeg:N]
  end
  # inds = @. inds_obs - ibeg + 1
  SM, SINK, Q
end


export cal_Q_Zeng2009!, RichardsEquation_Zeng2009
export solve_SM_Zeng2009

