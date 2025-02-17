# ! 注意
# - CoLM中，z向下为正
function soil_moisture_Zeng2009(soil::Soil{FT}, Q0::FT=0.0) where {FT<:Real}
  cal_θEψE!(soil)
  (; N, jwt, ibeg) = soil
  (; ψ, θ, K₊ₕ, ψE, sink) = soil
  (; θ_sat, param) = soil.param
  cal_K!(soil, θ)
  # cal_ψ!(soil, θ)

  # Q = cal_Q_Zeng2009!(soil, soil.θ)
  dt = soil.dt / 3600 # [s] -> [h]
  zwt = soil.zwt * 100 # [m] -> [cm]
  z = soil.z_cm
  Δz = soil.Δz_cm
  # Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100

  dKdθ = zeros(FT, N)
  dψdθ = zeros(FT, N + 1)
  ∂qᵢ∂θᵢ = zeros(FT, N + 1)
  ∂qᵢ∂θᵢ₊₁ = zeros(FT, N + 1)

  (; Q, a, b, c, d, e, f) = soil
  dθ = soil.du

  # Aquifer (11th) layer
  z[N+1] = 0.5 * (zwt + z[N]) # 动态调整最后一层的深度，高明! 中间位置
  Δz_gw = jwt >= N ? abs(zwt - z[N]) : Δz[N]

  # Hydraulic conductivity and soil matric potential and their derivatives
  for j in 1:N
    par = param[j]  # updated to ensure correct index usage

    j2 = min(N, j + 1)
    _θ = 0.5(θ[j] + θ[j2])
    _θsat = 0.5(θ_sat[j] + θ_sat[j2])
    se = clamp(_θ / _θsat, 0.01, 1.0)

    dKdθ[j] = Retention_∂K∂Se(se, par) / (2 * _θsat)  # CLM5, Eq. 7.87
    ψ[j] = Retention_ψ(θ[j], par)
    dψdθ[j] = Retention_∂ψ∂θ(ψ[j], par) # CLM5, Eq. 7.85
  end

  i = N
  if jwt == N
    se = 0.5 * (1.0 + θ[i] / θ_sat[i])
    se = clamp(se, 0.01, 1.0)
    # compute for aquifer layer [N+1]
    par = param[i]
    ψ[i+1] = Retention_ψ_Se(se, par) # N+1层的ψ，用的是第N层
    dψdθ[i+1] = Retention_∂ψ∂θ(ψ[i+1], par) #
  end

  for i in 1:N
    dz = (z[i+1] - z[i])
    dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
    Q[i] = -K₊ₕ[i] * dψ / dz

    ∂qᵢ∂θᵢ[i] = -(-K₊ₕ[i] * dψdθ[i] + dψ * dKdθ[i]) / dz
    ∂qᵢ∂θᵢ₊₁[i] = -(K₊ₕ[i] * dψdθ[i+1] + dψ * dKdθ[i]) / dz
  end
  
  i = N
  ∂qᵢ∂θᵢ[i+1] = 0.0
  ∂qᵢ∂θᵢ₊₁[i+1] = 0.0

  if jwt < N
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
  return Q
end


function error_SM(soil::Soil{FT}) where {FT<:Real}
  (; N, θ, θ_prev, Q, Q0) = soil
  dt = soil.dt / 3600 # [s] -> [h]
  Δz = soil.Δz_cm
  QN = Q[N]
  dθ = 0
  for i = 1:N
    dθ += (θ[i] - θ_prev[i]) * Δz[i] # in cm
  end
  obs = (QN - Q0) * dt # cm
  sim = dθ

  bias = sim - obs
  perc = bias / obs * 100 |> _round

  info = (; obs, sim, bias, perc, dθ, QN, Q0)
  info
end

_round(x::Real; digits=3) = round(x; digits)

export cal_θEψE!, soil_moisture_Zeng2009
