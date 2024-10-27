function soil_WaterFlux!(Q::V, θ::V, K::V, ψ::V, z::V;
  param::ParamVanGenuchten{T},
  ψ0::T=NaN, Q0::T=NaN, method="ψ0") where {T<:Real,V<:AbstractVector{T}}
  n = length(θ)
  @inbounds for i = 1:n
    K[i] = van_genuchten_K(θ[i]; param)
    ψ[i] = van_genuchten_ψ(θ[i]; param)
  end
  # @. K = van_genuchten_K(u; param)
  # @. ψ = van_genuchten_ψ(u; param)
  if method == "ψ0"
    Q0 = -K[1] * ((ψ0 - ψ[1]) / (0 - z[1]) + 1) # [cm/s]
  elseif method == "Q0"
    # Q0 = Q0
  end

  @inbounds for i in 1:n-1
    K₊ₕ = (K[i] + K[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    Q[i] = -K₊ₕ * ((ψ[i] - ψ[i+1]) / Δz₊ₕ + 1)
  end
  Q[n] = -K[n] # 尾部重力排水
  Q0
end


"""
# Arguments
- `method`: 
  + `ψ0`: ψ0 boundary condition, 第一类边界条件
  + `Q0`: Q0 boundary condition, 第二类边界条件
"""
function RichardsEquation(dθ::AbstractVector{T}, θ::AbstractVector{T}, p::Soil{T}, t; method="ψ0") where {T<:Real}
  p.timestep += 1
  (; ibeg, n, Q, K, ψ, ψ0, Q0, sink) = p # Δz, z, 
  Δz = p.Δz_cm
  z = p.z_cm

  Q0 = soil_WaterFlux!(Q, θ, K, ψ, z; param=p.param_water, ψ0, Q0, method)

  dθ[ibeg] = -((Q0 - Q[ibeg]) + sink[ibeg]) / Δz[ibeg]
  @inbounds for i in ibeg+1:n
    dθ[i] = -(Q[i-1] - Q[i]) / Δz[i] - sink[i] / Δz[i]
  end
end


function RichardsEquation_partial(dθ, θ, p::Soil, t; method="ψ0")
  (; ibeg) = p
  p.du[ibeg:end] .= dθ
  p.u[ibeg:end] .= θ
  RichardsEquation(p.du, p.u, p, t; method)
  dθ .= p.du[ibeg:end]
  return nothing
end
