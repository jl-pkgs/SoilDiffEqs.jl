# Calculate the equilibrium water content based on the water table depth
function cal_θEψE!(soil::Soil{T}) where {T<:Real}
  (; N, θE, ψE, z, zwt, jwt, method_retention) = soil
  (; param, θ_sat, ψ_sat) = soil.param
  iszero_ψsat = method_retention == "van_Genuchten" ? true : false

  # Δz = soil.Δz_cm
  z = soil.z_cm
  zwt = soil.zwt * 100
  soil.jwt = find_jwt(soil.z₊ₕ, soil.zwt) # ? 
  jwt = soil.jwt

  for j = 1:N
    z0 = z[j-1]
    z1 = z[j]
    par = param[j]
    _ψsat = iszero_ψsat ? 0.0 : ψ_sat[j]
    θE[j] = cal_θE(z1, z0, zwt, _ψsat, par)
    ψE[j] = Retention_ψ(θE[j], par)
  end

  # If zwt below soil column, calculate ψE for the 11th layer
  j = N
  if jwt == N
    # 积分在：z₊ₕ[N] ~ zwt，最后一层的θE、ψE代表的区间
    par = param[j]
    _ψsat = iszero_ψsat ? 0.0 : ψ_sat[j]
    θE[j+1] = cal_θE(zwt, z[j], zwt, _ψsat, par)
    ψE[j+1] = Retention_ψ(θE[j+1], par)
  elseif jwt < N # 最后一层饱和
    θE[j+1] = θ_sat[j]
    ψE[j+1] = T(0.0)
  end
  return ψE
end

# 中间变量
# ! 注意
# - CoLM中，z向下为正
function soil_moisture_zeng2009(soil::Soil{FT}, qflx_infl::FT=0.0) where {FT<:Real}
  (; N, zwt, θ, ψE) = soil
  (; θ_sat, param) = soil.param
  dt = soil.dt / 3600 # [s] -> [h]

  dKdθ = zeros(FT, N)
  dψdθ = zeros(FT, N)
  dqidθ0 = zeros(FT, N)
  dqidθ1 = zeros(FT, N)
  dqodθ1 = zeros(FT, N)
  dqodθ2 = zeros(FT, N)

  z = soil.z_cm
  Δz = soil.Δz_cm
  # Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100
  cal_θEψE!(soil)

  # Hydraulic conductivity and soil matric potential and their derivatives
  for j in 1:N
    par = param[j]  # updated to ensure correct index usage
    
    j2 = min(N, j + 1)
    _θ = 0.5(θ[j] + θ[j2])
    _θsat = 0.5(θ_sat[j] + θ_sat[j2])
    se = clamp(_θ / _θsat, 0.01, 1.0)

    # s2 = Ksat[j] * se^(2.0 * B[j] + 2.0)
    # dKdθ[j] = (2 * B[j] + 3.0) * s2 / (2 * _θsat)    
    dKdθ[j] = Retention_∂K∂Se(se, par) / (2 * _θsat)  # CLM5, Eq. 7.87
    ψ[j] = Retention_ψ(θ[j], par)
    dψdθ[j] = Retention_∂ψ∂θ(ψ[j], par) # CLM5, Eq. 7.85
  end

  # Set up r, a, b, and c vectors for tridiagonal solution
  i = 1
  qin[i] = qflx_infl # Infiltration rate, [mm H₂O/s]
  dz = (z[i+1] - z[i])
  dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
  qout[i] = -K₊ₕ[i] * dψ / dz
  dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
  dqodθ2[i] = -(K₊ₕ[i] * dψdθ[i+1] + dψ * dKdθ[i]) / dz

  sdamp = 0.0 # extrapolates θ dependence of evaporation
  rmx[i] = qin[i] - qout[i] - ET[i]
  amx[i] = 0.0
  bmx[i] = Δz[i] * (sdamp + 1.0 / dt) + dqodθ1[i] # ! -Δz[i]
  cmx[i] = dqodθ2[i]

  for i in 2:N-1
    dz = (z[i] - z[i-1]) # pos
    dψ = ψ[i] - ψE[i] - (ψ[i-1] - ψE[i-1])
    qin[i] = -K₊ₕ[i-1] * dψ / dz
    dqidθ0[i] = -(-K₊ₕ[i-1] * dψdθ[i-1] + dψ * dKdθ[i-1]) / dz
    dqidθ1[i] = -(K₊ₕ[i-1] * dψdθ[i] + dψ * dKdθ[i-1]) / dz

    dz = (z[i+1] - z[i])
    dψ = ψ[i+1] - ψE[i+1] - (ψ[i] - ψE[i])
    qout[i] = -K₊ₕ[i] * dψ / dz
    dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    dqodθ2[i] = -(K₊ₕ[i] * dψdθ[i+1] + dψ * dKdθ[i]) / dz

    amx[i] = -dqidθ0[i]
    bmx[i] = Δz[i] / dt - dqidθ1[i] + dqodθ1[i]
    cmx[i] = dqodθ2[i]
    rmx[i] = qin[i] - qout[i] - ET[i]
  end

  # Aquifer (11th) layer
  z_gw = 0.5 * (1e2 * zwt + z[N]) # 动态调整最后一层的深度，高明! 中间位置
  Δz_gw = Δz[N] # if jwt < N
  jwt >= N && (Δz_gw = 1e2 * zwt - z[N]) # in cm
  
  # Node j=N (bottom)
  i = N
  if jwt < N  # water table is in soil column
    dz = (z[i] - z[i-1])
    dψ = (ψ[i] - ψE[i]) - (ψ[i-1] - ψE[i-1])
    qin[i] = -K₊ₕ[i-1] * dψ / dz
    dqidθ0[i] = -(-K₊ₕ[i-1] * dψdθ[i-1] + dψ * dKdθ[i-1]) / dz
    dqidθ1[i] = -(K₊ₕ[i-1] * dψdθ[i] + dψ * dKdθ[i-1]) / dz
    qout[i] = 0.0
    dqodθ1[i] = 0.0

    amx[i] = -dqidθ0[i]
    bmx[i] = Δz[i] / dt - dqidθ1[i] + dqodθ1[i]
    cmx[i] = 0.0
    rmx[i] = qin[i] - qout[i] - ET[i]

    # next set up aquifer layer; hydrologically inactive
    amx[i+1] = 0.0
    bmx[i+1] = Δz_gw / dt # Δz[i+1] / dt
    cmx[i+1] = 0.0
    rmx[i+1] = 0.0
  else  # water table is below soil column
    # compute aquifer soil moisture as average of layer 10 and saturation
    ## N层
    se = 0.5 * (1.0 + θ[i] / θ_sat[i])
    se = clamp(se, 0.01, 1.0)

    # compute for aquifer layer [N+1]
    par = param[i]
    _ψ = Retention_ψ_se(se, par) # N+1层的ψ，用的是第N层
    dψdθ1 = Retention_∂ψ∂θ(_ψ, par) #
    # dψdθ1 = -B[i] * _ψ / (se * θ_sat[i])

    # first set up bottom layer of soil column
    dz = (z[i] - z[i-1])
    dψ = (ψ[i] - ψE[i]) - (ψ[i-1] - ψE[i-1])
    qin[i] = -K₊ₕ[i-1] * dψ / dz
    dqidθ0[i] = -(-K₊ₕ[i-1] * dψdθ[i-1] + dψ * dKdθ[i-1]) / dz
    dqidθ1[i] = -(K₊ₕ[i-1] * dψdθ[i] + dψ * dKdθ[i-1]) / dz

    # dz = (z[i+1] - z[i])
    dz = (z_gw - z[i])
    dψ = (_ψ - ψE[i+1]) - (ψ[i] - ψE[i])
    qout[i] = -K₊ₕ[i] * dψ / dz
    dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    dqodθ2[i] = -(K₊ₕ[i] * dψdθ1 + dψ * dKdθ[i]) / dz

    amx[i] = -dqidθ0[i]
    bmx[i] = Δz[i] / dt - dqidθ1[i] + dqodθ1[i]
    cmx[i] = dqodθ2[i]
    rmx[i] = qin[i] - qout[i] - ET[i]

    ## N+1层
    # next set up aquifer layer; dz/num unchanged, qin=qout
    qin[i+1] = qout[i]
    dqidθ0[i+1] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    dqidθ1[i+1] = -(K₊ₕ[i] * dψdθ1 + dψ * dKdθ[i]) / dz
    qout[i+1] = 0.0  # zero-flow bottom boundary condition
    dqodθ1[i+1] = 0.0  # zero-flow bottom boundary condition

    amx[i+1] = -dqidθ0[i+1]
    # bmx[i+1] = Δz[i+1] / dt - dqidθ1[i+1] + dqodθ1[i+1]
    bmx[i+1] = Δz_gw / dt - dqidθ1[i+1] + dqodθ1[i+1]
    cmx[i+1] = 0.0
    rmx[i+1] = qin[i+1] - qout[i+1]
  end
  Tridiagonal(amx, bmx, cmx, rmx, dθ; ibeg) # Solve for dθ

  # Renew the mass of liquid water also compute qcharge from dθ in aquifer layer
  # update in drainage for case jwt < N
  for j in 1:N
    θ[j] += dθ[j] * Δz[j]
  end
end

export cal_θEψE!
