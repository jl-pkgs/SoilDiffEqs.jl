"""
# Arguments
- `method`: 
  + `ψ0`: ψ0 boundary condition, 第一类边界条件
  + `Q0`: Q0 boundary condition, 第二类边界条件
"""
function RichardsEquation(dθ::AbstractVector{T}, u::AbstractVector{T}, p::Soil{T}, t; method="ψ0") where {T<:Real}
  p.timestep += 1
  (; n, Δz, z, Q, K, ψ, ψ0, sink) = p
  param = p.param_water
  @inbounds for i = 1:n
    K[i] = van_genuchten_K(u[i]; param)
    ψ[i] = van_genuchten_ψ(u[i]; param)
  end
  # @. K = van_genuchten_K(u; param)
  # @. ψ = van_genuchten_ψ(u; param)
  if method == "ψ0"
    Q0 = -K[1] * ((ψ0 - ψ[1]) / (0 - z[1]) + 1)
  elseif method == "Q0"
    Q0 = p.Q0
  end

  @inbounds for i in 1:n-1
    K₊ₕ = (K[i] + K[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    Q[i] = -K₊ₕ * ((ψ[i] - ψ[i+1]) / Δz₊ₕ + 1)
  end
  Q[n] = -K[n]

  dθ[1] = -(Q0 - Q[1]) / Δz[1] - sink[1] / Δz[1]
  @inbounds for i in 2:n
    dθ[i] = -(Q[i-1] - Q[i]) / Δz[i] - sink[i] / Δz[i]
  end
end
