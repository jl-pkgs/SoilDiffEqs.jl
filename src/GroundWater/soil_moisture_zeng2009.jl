function soil_moisture_zeng2009(soil::Soil{FT}) where {FT<:Real}
  dtime = get_step_size()

  # Compute jwt index
  jwt = find_jwt(zi, zwt)
  vwc_zwt = θ_sat[N]

  # ! 核心参考部分
  # Calculate the equilibrium water content based on the water table depth
  ψE = cal_θeψe!(θE, ψE, z, zwt, jwt; θ_sat, ψ_sat, bsw, ψmin)

  # Hydraulic conductivity and soil matric potential and their derivatives
  sdamp = 0.0
  for j in 1:N
    if origflag == 1
      se = (θ[j] + θ[min(N, j + 1)]) / ((θ_sat[j] + θ_sat[min(N, j + 1)]))
    else
      se = (vwc_liq[j] + vwc_liq[min(N, j + 1)]) / ((θ_sat[j] + θ_sat[min(N, j + 1)]))
    end
    se = min(1.0, se)

    # K₊ₕ[j] = imped[j] * se * s2
    dKdθ[j] = imped[j] * (2.0 * bsw[j] + 3.0) * s2 * (1.0 / (θ_sat[j] + θ_sat[min(N, j + 1)]))

    _θ = origflag == 1 ? θ[j] : vwc_liq[j]
    ψ[j] = cal_ψ(_θ[j], θ_sat[j], ψ_sat[j], bsw[j]; ψmin)
    dψdθ[j] = -bsw[j] * ψ[j] / (_θ)
  end

  # Aquifer (11th) layer
  zmm[N+1] = 0.5 * (1e3 * zwt + zmm[N]) # 动态调整最后一层的深度，高明!
  if jwt < N
    dzmm[N+1] = dzmm[N] # 这里需要核对
  else
    dzmm[N+1] = (1e3 * zwt - zmm[N])
  end

  # Set up r, a, b, and c vectors for tridiagonal solution
  i = 1
  qin[i] = qflx_infl
  dz = (zmm[i+1] - zmm[i])
  dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
  qout[i] = -K₊ₕ[i] * dψ / dz
  dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
  dqodθ2[i] = -(K₊ₕ[i] * dψdθ[i+1] + dψ * dKdθ[i]) / dz

  rmx[i] = qin[i] - qout[i] - ET[i]
  amx[i] = 0.0
  bmx[i] = dzmm[i] * (sdamp + 1.0 / dtime) + dqodθ1[i]
  cmx[i] = dqodθ2[i]

  # Nodes j=2 to j=N-1
  for i in 2:N-1
    dz = (zmm[i] - zmm[i-1])
    dψ = ψ[i] - ψE[i] - (ψ[i-1] - ψE[i-1])
    qin[i] = -K₊ₕ[i-1] * dψ / dz
    dqidθ0[i] = -(-K₊ₕ[i-1] * dψdθ[i-1] + dψ * dKdθ[i-1]) / dz
    dqidθ1[i] = -(K₊ₕ[i-1] * dψdθ[i] + dψ * dKdθ[i-1]) / dz

    dz = (zmm[i+1] - zmm[i])
    dψ = ψ[i+1] - ψE[i+1] - (ψ[i] - ψE[i])
    qout[i] = -K₊ₕ[i] * dψ / dz
    dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    dqodθ2[i] = -(K₊ₕ[i] * dψdθ[i+1] + dψ * dKdθ[i]) / dz

    amx[i] = -dqidθ0[i]
    bmx[i] = dzmm[i] / dtime - dqidθ1[i] + dqodθ1[i]
    cmx[i] = dqodθ2[i]
    rmx[i] = qin[i] - qout[i] - ET[i]
  end

  # Node j=N (bottom)
  i = N
  if jwt < N  # water table is in soil column
    dz = (zmm[i] - zmm[i-1])
    dψ = (ψ[i] - ψE[i]) - (ψ[i-1] - ψE[i-1])
    qin[i] = -K₊ₕ[i-1] * dψ / dz
    dqidθ0[i] = -(-K₊ₕ[i-1] * dψdθ[i-1] + dψ * dKdθ[i-1]) / dz
    dqidθ1[i] = -(K₊ₕ[i-1] * dψdθ[i] + dψ * dKdθ[i-1]) / dz
    qout[i] = 0.0
    dqodθ1[i] = 0.0

    amx[i] = -dqidθ0[i]
    bmx[i] = dzmm[i] / dtime - dqidθ1[i] + dqodθ1[i]
    cmx[i] = 0.0
    rmx[i] = qin[i] - qout[i] - ET[i]

    # next set up aquifer layer; hydrologically inactive
    # TODO: 漂亮！
    amx[i+1] = 0.0
    bmx[i+1] = dzmm[i+1] / dtime
    cmx[i+1] = 0.0
    rmx[i+1] = 0.0
  else  # water table is below soil column
    # compute aquifer soil moisture as average of layer 10 and saturation
    if origflag == 1
      se = 0.5 * (1.0 + θ[i] / θ_sat[i])
    else
      se = 0.5 * ((vwc_zwt + vwc_liq[i]) / θ_sat[i])
    end
    se = clamp(se, 0.01, 1.0)

    # compute for aquifer layer
    _ψ = cal_ψ(se, ψ_sat[i], bsw[i]; ψmin) # N+1层的ψ
    dψdθ1 = -bsw[i] * _ψ / (se * θ_sat[i])

    # first set up bottom layer of soil column
    dz = (zmm[i] - zmm[i-1])
    dψ = (ψ[i] - ψE[i]) - (ψ[i-1] - ψE[i-1])
    qin[i] = -K₊ₕ[i-1] * dψ / dz
    dqidθ0[i] = -(-K₊ₕ[i-1] * dψdθ[i-1] + dψ * dKdθ[i-1]) / dz
    dqidθ1[i] = -(K₊ₕ[i-1] * dψdθ[i] + dψ * dKdθ[i-1]) / dz

    dz = (zmm[i+1] - zmm[i])
    dψ = (_ψ - ψE[i+1]) - (ψ[i] - ψE[i])
    qout[i] = -K₊ₕ[i] * dψ / dz
    dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    dqodθ2[i] = -(K₊ₕ[i] * dψdθ1 + dψ * dKdθ[i]) / dz

    amx[i] = -dqidθ0[i]
    bmx[i] = dzmm[i] / dtime - dqidθ1[i] + dqodθ1[i]
    cmx[i] = dqodθ2[i]
    rmx[i] = qin[i] - qout[i] - ET[i]

    # next set up aquifer layer; dz/num unchanged, qin=qout
    qin[i+1] = qout[i]
    dqidθ0[i+1] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    dqidθ1[i+1] = -(K₊ₕ[i] * dψdθ1 + dψ * dKdθ[i]) / dz
    qout[i+1] = 0.0  # zero-flow bottom boundary condition
    dqodθ1[i+1] = 0.0  # zero-flow bottom boundary condition

    amx[i+1] = -dqidθ0[i+1]
    bmx[i+1] = dzmm[i+1] / dtime - dqidθ1[i+1] + dqodθ1[i+1]
    cmx[i+1] = 0.0
    rmx[i+1] = qin[i+1] - qout[i+1]
    # TODO: 漂亮！，最后一层不考虑蒸发。
  end
  Tridiagonal(amx, bmx, cmx, rmx, dθ; jtop=1) # Solve for dθ

  # Renew the mass of liquid water also compute qcharge from dθ in aquifer layer
  # update in drainage for case jwt < N
  for j in 1:N
    θ_liq[j] += dθ[j] * dzmm[j]
  end
end


# # calculate qcharge for case jwt < N
# if jwt < N
#   j0 = min(1, jwt)
#   j1 = jwt + 1

#   wh_zwt = 0.0
#   # Recharge rate qcharge to groundwater (positive to aquifer)

#   se = clamp(θ[jwt+1] / θ_sat[jwt+1], 0.01, 1.0)
#   # scs: this is the expression for unsaturated K
#   _K = imped[jwt+1] * Ksat[jwt+1] * se^(2.0 * bsw[jwt+1] + 3.0)

#   _ψ = max(ψmin, ψ[j0])
#   wh = _ψ - ψE[j0]  # 这里是向地下水的排泄，Zeng2009, Eq.14

#   if jwt == 0
#     qcharge = -_K * (wh_zwt - wh) / ((zwt + 1.0e-3) * 1000.0)
#   else
#     qcharge = -_K * (wh_zwt - wh) / ((zwt - z[jwt]) * 1000.0 * 2.0)
#   end
#   # To limit qcharge (for the first several timesteps)
#   qcharge = max(qcharge, -10.0 / dtime, 10.0 / dtime)
# else
#   # if water table is below soil column, compute qcharge from dθ(11)
#   qcharge = dθ[N+1] * dzmm[N+1] / dtime
# end

# # compute the water deficit and reset negative liquid water content
# qflx_deficit = 0.0
# for j in 1:N
#   if θ_liq[j] < 0.0
#     qflx_deficit -= θ_liq[j]
#   end
# end
