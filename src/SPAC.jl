## 主要模拟器

"""
    solve_SM_ODE(soil, Tsurf; solver)

solver = Tsit5()
solver = Rosenbrock23()
solver = Rodas5(autodiff=false)  
"""
function solve_SM_Zeng2009(soil; solver, verbose=false,
  reltol=1e-3, abstol=1e-3, ϵ=1e-4,
  GW=true,
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
  Z = zeros(ntime)

  SM[1, :] .= soil.θ[ibeg:N]

  nchunk = 50
  chunksize = ceil(Int, ntime / nchunk)
  nchunk = ceil(Int, ntime / chunksize)
  p = Progress(nchunk)

  zwt0 = Z[1] = soil.zwt
  # 加入土壤含水量限制因子
  for i = 2:ntime
    (verbose && mod(i, chunksize) == 0) && next!(p)
    # soil.θ0 = θ_surf[i]
    prob.u0 .= soil.θ[ibeg:N]

    k = 2
    β = soil.θ[k] / soil.param.θ_sat[k]
    ## 这里要添加根系抽水过程
    soil.sink[k] = ET[i] / 10 * β # mm/h to cm/h，目前根系抽水写在了第二层

    soil.θ_prev[1:N] .= soil.θ[1:N]
    sol = _solve(prob, solver; reltol, abstol, saveat=dt)
    soil.θ[ibeg:N] .= sol.u[end] # 更新这个时刻的结果

    ## 之后调整地下水的水位
    ∑ = -soil.Q[N] * 10 # [cm h-1] to [mm]
    soil.zwt, soil.wa, soil.uex = GW_Update_ZWT!(soil, soil.θ; ∑)
    uex, ∑ = GW_Correctθ!(soil, soil.θ; exceed2surf=false)
    if abs(∑) >= ϵ
      ## 再调整一次
      soil.zwt, soil.wa, soil.uex = GW_Update_ZWT!(soil, soil.θ; ∑)
    end

    !GW && (soil.zwt = zwt0) # 关闭地下水模块
    Z[i] = soil.zwt

    SINK[i, :] .= soil.sink[1:N]
    Q[i, :] = cal_Q_Zeng2009!(soil, soil.θ)
    SM[i, :] .= soil.θ[ibeg:N]
  end
  # inds = @. inds_obs - ibeg + 1
  SM, SINK, Q, Z
end

export solve_SM_Zeng2009
