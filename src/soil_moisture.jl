import HydroTools: soil_moisture!

# soil_moisture!(θ, ψ, sink, ψ0, dz, dt, param)
function soil_moisture!(θ::V, ψ::V, sink::V, ψ0::T, dz::V, dt, param; fun=van_Genuchten) where {
  T<:Real,V<:AbstractVector{T}}
  z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  n = length(dz)

  θ_n = deepcopy(θ)
  ψ_n = deepcopy(ψ)
  isnothing(sink) && (sink = zeros(n))

  # θ = zeros(n)
  K = zeros(n)
  Cap = zeros(n)
  ψ_pred = zeros(n)
  K₊ₕ = zeros(n)

  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ[i]; param)
  end

  for i = 1:n-1
    K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end

  K0₊ₕ = K[1]
  dz0₊ₕ = 0.5 * dz[1]

  a = zeros(n)
  c = zeros(n)
  b = zeros(n)
  d = zeros(n)

  @inbounds for i = 1:n
    if i == 1
      a[i] = 0
      c[i] = -K₊ₕ[i] / dz₊ₕ[i]
      b[i] = Cap[i] * dz[i] / (0.5 * dt) + K0₊ₕ / dz0₊ₕ - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K0₊ₕ / dz0₊ₕ * ψ0 + K0₊ₕ - K₊ₕ[i]
    elseif i < n
      a[i] = -K₊ₕ[i-1] / dz₊ₕ[i-1]
      c[i] = -K₊ₕ[i] / dz₊ₕ[i]
      b[i] = Cap[i] * dz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[i-1] - K₊ₕ[i]
    elseif i == n
      a[i] = -K₊ₕ[n-1] / dz₊ₕ[n-1]
      c[i] = 0
      b[i] = Cap[i] * dz[i] / (0.5 * dt) - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / (0.5 * dt) * ψ[i] + K₊ₕ[n-1] - K[i]
    end
    d[i] -= sink[i]
  end
  ψ_pred .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1/2 time

  ## update: θ, K and Cap
  for i in 1:n
    θ[i], K[i], Cap[i] = fun(ψ_pred[i]; param)
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
      c[i] = -K₊ₕ[i] / (2 * dz₊ₕ[i])
      b[i] = Cap[i] * dz[i] / dt - c[i] + K0₊ₕ / (2 * dz0₊ₕ)
      d[i] = Cap[i] * dz[i] / dt * ψ[i] +
             K0₊ₕ / (2 * dz0₊ₕ) * (2ψ0 - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K0₊ₕ - K₊ₕ[i]
    elseif i < n
      a[i] = -K₊ₕ[i-1] / (2 * dz₊ₕ[i-1])
      c[i] = -K₊ₕ[i] / (2 * dz₊ₕ[i])
      b[i] = Cap[i] * dz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) +
             c[i] * (ψ[i] - ψ[i+1]) + K₊ₕ[i-1] - K₊ₕ[i]
    else
      i == n
      a[i] = -K₊ₕ[i-1] / (2 * dz₊ₕ[i-1])
      c[i] = 0
      b[i] = Cap[i] * dz[i] / dt - a[i] - c[i]
      d[i] = Cap[i] * dz[i] / dt * ψ[i] - a[i] * (ψ[i-1] - ψ[i]) + K₊ₕ[i-1] - K[i]
    end
    d[i] -= sink[i]
  end
  ψ .= tridiagonal_solver(a, b, c, d) # Solve for ψ at n+1
  # for i in 1:n
  #   θ[i], K[i], Cap[i] = fun(ψ_pred[i]; param)
  # end
  
  # --- Check water balance
  Q0 = -K0₊ₕ / (2 * dz0₊ₕ) * ((ψ0 - ψ_n[1]) + (ψ0 - ψ[1])) - K0₊ₕ # 两个时刻的
  QN = -K[n]

  dθ = 0
  for i = 1:n
    dθ += (θ[i] - θ_n[i]) * dz[i]
  end

  err = dθ - (QN - Q0) * dt
  Q0, QN, dθ, err
end
