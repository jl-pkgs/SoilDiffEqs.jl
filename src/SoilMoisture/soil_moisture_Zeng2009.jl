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
  jwt < N && (Q[N] = 0.0)

  i = 1
  sdamp = 0.0 # extrapolates θ dependence of evaporation
  a[i] = 0.0
  b[i] = Δz[i] * (sdamp - 1.0 / dt) + ∂qᵢ∂θᵢ[i] # 
  c[i] = ∂qᵢ∂θᵢ₊₁[i]
  d[i] = Q0 - Q[i] + sink[i]

  for i in 2:N-1
    a[i] = -∂qᵢ∂θᵢ[i-1]
    b[i] = -Δz[i] / dt - ∂qᵢ∂θᵢ₊₁[i-1] + ∂qᵢ∂θᵢ[i]
    c[i] = ∂qᵢ∂θᵢ₊₁[i]
    d[i] = Q[i-1] - Q[i] + sink[i]
  end

  # Node j=N (bottom)
  i = N
  if jwt < N  # water table is in soil column
    # qi = 0.0, ∂qᵢ∂θᵢ[i] = 0.0
    a[i] = -∂qᵢ∂θᵢ[i-1]
    b[i] = -Δz[i] / dt - ∂qᵢ∂θᵢ₊₁[i-1]
    c[i] = 0.0
    d[i] = Q[i-1] - Q[i] + sink[i]

    # next set up aquifer layer; hydrologically inactive
    a[i+1] = 0.0
    b[i+1] = Δz_gw / dt # Δz[i+1] / dt
    c[i+1] = 0.0
    d[i+1] = 0.0
  else  # water table is below soil column
    # compute aquifer soil moisture as average of layer 10 and saturation
    a[i] = -∂qᵢ∂θᵢ[i-1]
    b[i] = -Δz[i] / dt - ∂qᵢ∂θᵢ₊₁[i-1] + ∂qᵢ∂θᵢ[i]
    c[i] = ∂qᵢ∂θᵢ₊₁[i]
    d[i] = Q[i-1] - Q[i] + sink[i]
    ## N+1层
    # next set up aquifer layer; dz/num unchanged, qin=Q
    # q[i+1] = 0.0  # zero-flow bottom boundary condition
    a[i+1] = -∂qᵢ∂θᵢ[i]
    b[i+1] = -Δz_gw / dt - ∂qᵢ∂θᵢ₊₁[i]
    c[i+1] = 0.0
    d[i+1] = Q[i] # 最后一层，地下水不考虑蒸发
  end
  tridiagonal_solver!(a, b, c, d, e, f, dθ; ibeg, N=N + 1)

  # Renew the mass of liquid water also compute qcharge from dθ in aquifer layer
  # update in drainage for case jwt < N
  for j in 1:N
    θ[j] += dθ[j] * Δz[j]
  end
  return Q
end


export cal_θEψE!, soil_moisture_Zeng2009
