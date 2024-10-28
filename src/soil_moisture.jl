function update_K₊ₕ!(soil)
  (; n, ibeg, z, z₊ₕ, K, K₊ₕ) = soil
  for i = ibeg:n-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    K₊ₕ[i] = K[i] * K[i+1] * (d1 + d2) / (K[i] * d2 + K[i+1] * d1) # Eq. 5.16, 
    # K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
end

# soil_moisture!(soil, sink, ψ0, param)
function soil_moisture!(soil::Soil, sink::V, ψ0::T;
  fun=van_Genuchten) where {T<:Real,V<:AbstractVector{T}}

  isnothing(sink) && (sink = zeros(n))

  (; n, dt, #Δz, Δz₊ₕ,
    ψ, ibeg,
    θ, Cap, K, ψ_next, K₊ₕ, θ_prev, ψ_prev, a, b, c, d) = soil
  Δz = soil.Δz_cm
  Δz₊ₕ = soil.Δz₊ₕ_cm
  param = soil.param_water

  θ_prev .= θ # backup
  ψ_prev .= ψ

  for i in ibeg:n
    θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  end
  update_K₊ₕ!(soil)

  K0₊ₕ = K[ibeg]
  dz0₊ₕ = ibeg == 1 ? 0.5 * Δz[1] : Δz₊ₕ[ibeg-1]
  # dz0₊ₕ = 0.5 * Δz[1] # ? 

  @inbounds for i = ibeg:n
    if i == ibeg
      a[i] = 0
      c[i] = -K₊ₕ[i] / Δz₊ₕ[i]
      b[i] = Cap[i] * Δz[i] / (0.5 * dt) + K0₊ₕ / dz0₊ₕ - c[i]
      d[i] = Cap[i] * Δz[i] / (0.5 * dt) * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
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

  _a = @view a[ibeg:end]
  _b = @view b[ibeg:end]
  _c = @view c[ibeg:end]
  _d = @view d[ibeg:end]
  ψ_next[ibeg:end] .= tridiagonal_solver(_a, _b, _c, _d)
  # ψ_next .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1/2 time

  ## update: θ, K and Cap
  for i in ibeg:n
    θ[i], K[i], Cap[i] = fun(ψ_next[i]; param)
  end
  update_K₊ₕ!(soil)
  K0₊ₕ = K[ibeg] # 可以按照同样的方法，设置

  ## second round
  # Terms for tridiagonal matrix
  @inbounds for i = ibeg:n
    if i == ibeg
      a[i] = 0
      c[i] = -K₊ₕ[i] / (2 * Δz₊ₕ[i])
      b[i] = Cap[i] * Δz[i] / dt - c[i] + K0₊ₕ / (2 * dz0₊ₕ)
      d[i] = Cap[i] * Δz[i] / dt * ψ[i] +
             K0₊ₕ / (2 * dz0₊ₕ) * (2ψ0 - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K0₊ₕ - K₊ₕ[i]
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
  _a = @view a[ibeg:end]
  _b = @view b[ibeg:end]
  _c = @view c[ibeg:end]
  _d = @view d[ibeg:end]
  ψ[ibeg:end] .= tridiagonal_solver(_a, _b, _c, _d)
  # ψ .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1

  ## variables already updated by dot operation
  # for i in ibeg:n
  #   θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  # end

  ## Check water balance
  Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_prev[1]) + (ψ0 - ψ[1])) - K0₊ₕ # 两个时刻的
  QN = -K[n]

  dθ = 0
  for i = 1:n
    dθ += (θ[i] - θ_prev[i]) * Δz[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
