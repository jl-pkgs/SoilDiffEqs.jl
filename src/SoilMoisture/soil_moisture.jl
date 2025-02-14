# soil_moisture!(soil, sink, ψ0, param)
function soil_moisture!(soil::Soil, sink::V, ψ0::T;
  debug::Bool=true) where {T<:Real,V<:AbstractVector{T}}

  (; N, #Δz, Δz₊ₕ,
    ψ, ibeg,
    θ, ∂θ∂ψ, K, ψ_next, K₊ₕ, θ_prev, ψ_prev,
    a, b, c, d, e, f) = soil
  dt = soil.dt / 3600 # [s] -> [h]
  Δz = soil.Δz_cm
  Δz₊ₕ = soil.Δz₊ₕ_cm

  θ_prev .= θ # backup
  ψ_prev .= ψ

  # cal_ψ!(soil, θ)     # ψ, 以θ为导向
  cal_θ!(soil, ψ)
  cal_K!(soil, θ)     # K
  cal_∂θ∂ψ!(soil, ψ)  # ∂θ∂ψ
  
  i0 = max(ibeg - 1, 1)
  K0₊ₕ = ibeg == 1 ? K[1] : K₊ₕ[i0]
  dz0₊ₕ = ibeg == 1 ? 0.5 * Δz[1] : Δz₊ₕ[i0]
  # dz0₊ₕ = 0.5 * Δz[1] # ? 
  dt_half = 0.5 * dt

  # first round:
  @inbounds for i = ibeg:N
    if i == ibeg
      a[i] = 0.0
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      b[i] = ∂θ∂ψ[i] * Δz[i] / dt_half + K0₊ₕ / dz0₊ₕ - c[i]
      d[i] = ∂θ∂ψ[i] * Δz[i] / dt_half * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
    elseif i < N
      a[i] = -K₊ₕ[i-1] / Δz₊ₕ[i-1]
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      b[i] = ∂θ∂ψ[i] * Δz[i] / dt_half - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz[i] / dt_half * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == N
      a[i] = -K₊ₕ[N-1] / Δz₊ₕ[N-1]
      c[i] = 0.0
      b[i] = ∂θ∂ψ[i] * Δz[i] / dt_half - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz[i] / dt_half * ψ[i] + K₊ₕ[N-1] - K[i]
    end
    d[i] -= sink[i]
  end
  tridiagonal_solver!(a, b, c, d, e, f, ψ_next; ibeg)

  ## update: θ, K and ∂θ∂ψ
  cal_θ!(soil, ψ_next)    # θ
  cal_K!(soil, θ)         # K
  cal_∂θ∂ψ!(soil, ψ_next) # ∂θ∂ψ
  K0₊ₕ = ibeg == 1 ? K[1] : K₊ₕ[i0]

  ## second round: in half step
  @inbounds for i = ibeg:N
    if i == ibeg
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ[i])
      b[i] = ∂θ∂ψ[i] * Δz[i] / dt - c[i] + K0₊ₕ / (2 * dz0₊ₕ)
      d[i] = ∂θ∂ψ[i] * Δz[i] / dt * ψ[i] +
             K0₊ₕ / (2 * dz0₊ₕ) * (2ψ0 - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K0₊ₕ - K₊ₕ[i]
    elseif i < N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ[i-1])
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ[i])
      b[i] = ∂θ∂ψ[i] * Δz[i] / dt - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ[i-1])
      c[i] = 0
      b[i] = ∂θ∂ψ[i] * Δz[i] / dt - a[i] - c[i]
      d[i] = ∂θ∂ψ[i] * Δz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
    d[i] -= sink[i]
  end
  tridiagonal_solver!(a, b, c, d, e, f, ψ; ibeg)
  # cal_θ!(soil, ψ)

  ## Check water balance
  Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_prev[1]) + (ψ0 - ψ[1])) - K0₊ₕ # 两个时刻的
  QN = -K[N]

  dθ = 0
  for i = 1:N
    dθ += (θ[i] - θ_prev[i]) * Δz[i]
  end

  err = dθ - (QN - Q0) * dt
  debug && return Q0, QN, dθ, err
  return nothing
end
