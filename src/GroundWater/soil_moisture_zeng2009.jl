# 中间变量

# ! 注意
# - CoLM中，z向下为正
function soil_moisture_zeng2009(soil::Soil{FT}, qflx_infl, dKdθ, dψdθ, dqidθ0, dqidθ1, dqodθ1, dqodθ2, dzmm) where {FT<:Real}
  (; N, dt, zwt, 
    z₊ₕ, z_cm, Δz_cm, Δz₊ₕ_cm) = soil
  
  dKdθ = zeros(FT, N)
  dψdθ = zeros(FT, N)
  dqidθ0 = zeros(FT, N)
  dqidθ1 = zeros(FT, N)
  dqodθ1 = zeros(FT, N)
  dqodθ2 = zeros(FT, N)
  
  # z_cm::Vector{FT} = z * 100         # cm, 向下为负
  # Δz_cm::Vector{FT} = Δz * 100
  # Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100
  jwt = find_jwt(z₊ₕ, zwt)

  ## 首先计算均衡水势
  # Calculate the equilibrium water content based on the water table depth
  ψE = cal_θeψe!(θE, ψE, z, zwt, jwt; θ_sat, ψ_sat, B, ψmin)

  # Hydraulic conductivity and soil matric potential and their derivatives
  for j in 1:N
    se = (θ[j] + θ[min(N, j + 1)]) / ((θ_sat[j] + θ_sat[min(N, j + 1)]))
    se = min(1.0, se)

    s2 = Ksat[j] * se^(2.0 * B[j] + 2.0)
    dKdθ[j] = (2B[j] + 3.0) * s2 * (1.0 / (θ_sat[j] + θ_sat[min(N, j + 1)]))    # CLM5, Eq. 7.87

    ψ[j] = cal_ψ(θ[j], θ_sat[j], ψ_sat[j], B[j]; ψmin)
    dψdθ[j] = -B[j] * ψ[j] / (_θ)                                               # CLM5, Eq. 7.85
  end

  sdamp = 0.0 # extrapolates θ dependence of evaporation

  # Set up r, a, b, and c vectors for tridiagonal solution
  i = 1
  qin[i] = qflx_infl # Infiltration rate, [mm H₂O/s]
  dz = (zmm[i+1] - zmm[i])
  dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
  qout[i] = -K₊ₕ[i] * dψ / dz
  dqodθ1[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
  dqodθ2[i] = -(K₊ₕ[i] * dψdθ[i+1] + dψ * dKdθ[i]) / dz

  rmx[i] = qin[i] - qout[i] - ET[i]
  amx[i] = 0.0
  bmx[i] = dzmm[i] * (sdamp + 1.0 / dt) + dqodθ1[i] # ! -dzmm[i]
  cmx[i] = dqodθ2[i]

  # Nodes j=2 to j=N-1
  for i in 2:N-1
    dz = (zmm[i] - zmm[i-1]) # pos
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
    bmx[i] = dzmm[i] / dt - dqidθ1[i] + dqodθ1[i]
    cmx[i] = dqodθ2[i]
    rmx[i] = qin[i] - qout[i] - ET[i]
  end

  # Aquifer (11th) layer
  zmm[N+1] = 0.5 * (1e3 * zwt + zmm[N]) # 动态调整最后一层的深度，高明! 中间位置
  if jwt < N
    dzmm[N+1] = dzmm[N] # 这里需要核对，没用到，无所谓
  else
    dzmm[N+1] = (1e3 * zwt - zmm[N])
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
    bmx[i] = dzmm[i] / dt - dqidθ1[i] + dqodθ1[i]
    cmx[i] = 0.0
    rmx[i] = qin[i] - qout[i] - ET[i]

    # next set up aquifer layer; hydrologically inactive
    amx[i+1] = 0.0
    bmx[i+1] = dzmm[i+1] / dt
    cmx[i+1] = 0.0
    rmx[i+1] = 0.0
  else  # water table is below soil column
    # compute aquifer soil moisture as average of layer 10 and saturation
    ## N层
    se = 0.5 * (1.0 + θ[i] / θ_sat[i])
    se = clamp(se, 0.01, 1.0)

    # compute for aquifer layer
    _ψ = cal_ψ(se, ψ_sat[i], B[i]; ψmin) # N+1层的ψ，用的是第N层
    dψdθ1 = -B[i] * _ψ / (se * θ_sat[i])

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
    bmx[i] = dzmm[i] / dt - dqidθ1[i] + dqodθ1[i]
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
    bmx[i+1] = dzmm[i+1] / dt - dqidθ1[i+1] + dqodθ1[i+1]
    cmx[i+1] = 0.0
    rmx[i+1] = qin[i+1] - qout[i+1]
  end
  Tridiagonal(amx, bmx, cmx, rmx, dθ; jtop=1) # Solve for dθ

  # Renew the mass of liquid water also compute qcharge from dθ in aquifer layer
  # update in drainage for case jwt < N
  for j in 1:N
    θ[j] += dθ[j] * dzmm[j]
  end
end
