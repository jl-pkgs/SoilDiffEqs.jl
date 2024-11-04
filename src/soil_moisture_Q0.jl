function soil_moisture_Q0!(soil::Soil{FT}, sink::V, Q0::FT;) where {
  FT<:Real,V<:AbstractVector{FT}}

  (; N, dt, #Δz, Δz₊ₕ,
    ψ,
    θ, ψ_next, Cap, K, K₊ₕ, θ_prev, ψ_prev,
    a, b, c, d, e, f) = soil
  Δz = soil.Δz_cm
  Δz₊ₕ = soil.Δz₊ₕ_cm

  θ_prev .= θ # backup
  ψ_prev .= ψ

  cal_θKCap!(soil, ψ)
  cal_K₊ₕ!(soil)

  K0₊ₕ = K[1]
  # dz0₊ₕ = 0.5 * Δz[1]
  @inbounds for i = 1:N
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      # b[i] = Cap[i] * Δz[i] / (0.5 * dt) + K0₊ₕ / dz0₊ₕ - c[i]
      # d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] - K₊ₕ[i] - Q0
    elseif i < N
      a[i] = -K₊ₕ[i-1] / Δz₊ₕ[i-1]
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == N
      a[i] = -K₊ₕ[N-1] / Δz₊ₕ[N-1]
      c[i] = 0
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[N-1] - K[i]
    end
    d[i] -= sink[i]
  end

  # Solve for ψ at N+1/2 time
  tridiagonal_solver!(a, b, c, d, e, f, ψ_next)
  # ψ_next .= tridiagonal_solver(a, b, c, d)

  ## update: θ, K and Cap
  cal_θKCap!(soil, ψ_next)
  cal_K₊ₕ!(soil)
  K0₊ₕ = K[1] # 可以按照同样的方法，设置

  ## second round
  # Terms for tridiagonal matrix
  @inbounds for i = 1:N
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ[i])
      # b[i] = Cap[i] * Δz[i] / dt + K0₊ₕ / (2 * dz0₊ₕ) - c[i]
      # d[i] = Cap[i] * Δz[i] / dt * ψ[i] + K0₊ₕ / (2 * dz0₊ₕ) * ψ0 +
      #        K0₊ₕ / (2 * dz0₊ₕ) * (ψ0 - ψ[i]) +
      #        c[i] * (ψ[i] - ψ[i+1]) + K0₊ₕ - K₊ₕ[i]
      b[i] = Cap[i] * Δz[i] / dt - c[i]
      d[i] = Cap[i] * Δz[i] / dt * ψ[i] - Q0 +
             c[i] * (ψ[i] - ψ[i+1]) - K₊ₕ[i]
    elseif i < N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ[i-1])
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ[i])
      b[i] = Cap[i] * Δz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    else
      i == N
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ[i-1])
      c[i] = 0
      b[i] = Cap[i] * Δz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
    d[i] -= sink[i]
  end
  # Solve for ψ at N+1
  tridiagonal_solver!(a, b, c, d, e, f, ψ)
  # ψ .= tridiagonal_solver(a, b, c, d)

  ## 获取ψ，之后更新θ
  # for i in 1:N
  #   θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  # end

  # --- Check water balance
  # Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_n[1]) + (ψ0 - ψ[1])) - K0₊ₕ
  QN = -K[N]
  dθ = 0
  for i = 1:N
    dθ += (θ[i] - θ_prev[i]) * Δz[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
