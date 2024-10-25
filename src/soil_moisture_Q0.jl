function soil_moisture_Q0!(soil::Soil{FT}, sink::V, Q0::FT, param; fun=van_Genuchten) where {
  FT<:Real,V<:AbstractVector{FT}}

  (; n, dt, Δz, Δz₊ₕ,
    θ, ψ, ψ_next, Cap, K, K₊ₕ, θ_prev, ψ_prev,
    a, b, c, d) = soil

  θ_prev .= θ # backup
  ψ_prev .= ψ

  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  end
  for i = 1:n-1
    K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end

  K0₊ₕ = K[1]
  # dz0₊ₕ = 0.5 * Δz[1]
  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      # b[i] = Cap[i] * Δz[i] / (0.5 * dt) + K0₊ₕ / dz0₊ₕ - c[i]
      # d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] - K₊ₕ[i] - Q0
    elseif i < n
      a[i] = -K₊ₕ[i-1] / Δz₊ₕ[i-1]
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == n
      a[i] = -K₊ₕ[n-1] / Δz₊ₕ[n-1]
      c[i] = 0
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[n-1] - K[i]
    end
    d[i] -= sink[i]
  end
  ψ_next .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1/2 time

  ## update: θ, K and Cap
  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ_next[i]; param)
  end
  for i = 1:n-1
    K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
  K0₊ₕ = K[1] # 可以按照同样的方法，设置

  ## second round
  # Terms for tridiagonal matrix
  @inbounds for i = 1:n
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
    elseif i < n
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ[i-1])
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ[i])
      b[i] = Cap[i] * Δz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    else
      i == n
      a[i] = -K₊ₕ[i-1] / (2 * Δz₊ₕ[i-1])
      c[i] = 0
      b[i] = Cap[i] * Δz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * Δz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
    d[i] -= sink[i]
  end
  ψ .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1
  ## 获取ψ，之后更新θ
  # for i in 1:n
  #   θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  # end

  # --- Check water balance
  # Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_n[1]) + (ψ0 - ψ[1])) - K0₊ₕ
  QN = -K[n]
  dθ = 0
  for i = 1:n
    dθ += (θ[i] - θ_prev[i]) * Δz[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
