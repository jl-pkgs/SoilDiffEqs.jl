# ! 注意
# - CoLM中，z向下为正
function soil_moisture_Zeng2009(soil::Soil{FT}, Q0::FT=0.0; ∂K₊ₕ∂θ=nothing) where {FT<:Real}
  cal_θEψE!(soil)
  (; N, jwt, ibeg) = soil
  (; ψ, θ, K₊ₕ, ψE, sink) = soil
  (; θ_sat, θ_res, param) = soil.param
  θ_res .= 0.0
  cal_K!(soil, θ)
  # cal_ψ!(soil, θ)

  # Q = cal_Q_Zeng2009!(soil, soil.θ)
  dt = soil.dt / 3600 # [s] -> [h]
  zwt = soil.zwt * 100 # [m] -> [cm]
  z = soil.z_cm
  Δz = soil.Δz_cm
  # Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100

  # ∂K₊ₕ∂θ = zeros(FT, N)
  ∂ψ∂θ = zeros(FT, N + 1)
  ∂qᵢ∂θᵢ = zeros(FT, N + 1)
  ∂qᵢ∂θᵢ₊₁ = zeros(FT, N + 1)

  (; Q, a, b, c, d, e, f) = soil
  dθ = soil.du

  # Aquifer (11th) layer
  z[N+1] = 0.5 * (zwt + z[N]) # 动态调整最后一层的深度，高明! 中间位置
  Δz_gw = jwt >= N ? abs(zwt - z[N]) : Δz[N]

  # Hydraulic conductivity and soil matric potential and their derivatives
  for i in 1:N
    par = param[i]  # updated to ensure correct index usage

    i2 = min(N, i + 1)
    _θ = 0.5(θ[i] + θ[i2])
    _θsat = 0.5(θ_sat[i] + θ_sat[i2])
    _θres = 0.5(θ_res[i] + θ_res[i2])
    se = clamp((_θ - _θres) / (_θsat - _θres), 0.01, 1.0)
    # se = clamp(_θ / _θsat, 0.01, 1.0)
    # ∂K₊ₕ∂θ[i] = Retention_∂K∂Se(se, par) / (2 * (_θsat - _θres))  # CLM5, Eq. 7.87
    ψ[i] = Retention_ψ(θ[i], par)
    ∂ψ∂θ[i] = Retention_∂ψ∂θ(ψ[i], par) # CLM5, Eq. 7.85
  end

  i = N
  if jwt == N
    _θ = 0.5 * (θ[i] + θ_sat[i])
    se = clamp((_θ - θ_res[i]) / (θ_sat[i] - θ_res[i]), 0.01, 1.0)
    # se = 0.5 * (1.0 + θ[i] / θ_sat[i])
    # compute for aquifer layer [N+1]
    par = param[i]
    ψ[i+1] = Retention_ψ_Se(se, par) # N+1层的ψ，用的是第N层
    ∂ψ∂θ[i+1] = Retention_∂ψ∂θ(ψ[i+1], par) #
  end

  # 如果告诉他正确的∂ψ∂θ, ∂K₊ₕ∂θ，它能否解对？
  for i in 1:N
    dz = (z[i+1] - z[i])
    dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
    Q[i] = -K₊ₕ[i] * dψ / dz

    ∂qᵢ∂θᵢ[i] = -(-K₊ₕ[i] * ∂ψ∂θ[i] + dψ * ∂K₊ₕ∂θ[i]) / dz
    ∂qᵢ∂θᵢ₊₁[i] = -(K₊ₕ[i] * ∂ψ∂θ[i+1] + dψ * ∂K₊ₕ∂θ[i]) / dz
  end
  
  i = N
  ∂qᵢ∂θᵢ[i+1] = 0.0
  ∂qᵢ∂θᵢ₊₁[i+1] = 0.0

  if jwt < N # 地下水侵入土壤
    Q[i] = 0.0
    ∂qᵢ∂θᵢ[i] = 0.0
    ∂qᵢ∂θᵢ₊₁[i] = 0.0
  end

  i = 1
  sdamp = 0.0 # extrapolates θ dependence of evaporation
  a[i] = 0.0
  b[i] = ∂qᵢ∂θᵢ[i] - Δz[i] / dt + Δz[i] * sdamp
  c[i] = ∂qᵢ∂θᵢ₊₁[i]
  d[i] = Q0 - Q[i] + sink[i]

  for i in 2:N
    a[i] = -∂qᵢ∂θᵢ[i-1]
    b[i] = ∂qᵢ∂θᵢ[i] - ∂qᵢ∂θᵢ₊₁[i-1] - Δz[i] / dt
    c[i] = ∂qᵢ∂θᵢ₊₁[i]
    d[i] = Q[i-1] - Q[i] + sink[i]
  end

  ## N+1 layer, qᵢ₊₁ = 0.0, qᵢ₊₁项移除
  i = N
  a[i+1] = -∂qᵢ∂θᵢ[i]
  b[i+1] = -∂qᵢ∂θᵢ₊₁[i] - Δz_gw / dt
  c[i+1] = 0.0
  d[i+1] = Q[i] # 最后一层，地下水不考虑蒸发

  tridiagonal_solver!(a, b, c, d, e, f, dθ; ibeg, N=N + 1)
  for j in 1:N
    θ[j] += dθ[j] * Δz[j] # update θ
  end
  return (; Q, ∂K₊ₕ∂θ, ∂ψ∂θ, ∂qᵢ∂θᵢ, ∂qᵢ∂θᵢ₊₁)
end


export cal_θEψE!, soil_moisture_Zeng2009
