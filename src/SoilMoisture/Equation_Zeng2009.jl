"""
- Q0: Infiltration rate, [cm h-1], should be negative or zero
"""
function cal_Q_Zeng2009!(soil::Soil{T}, θ::AbstractVector{T}; Q0::T=NaN) where {T<:Real}
  (; N, Q, ψ, ψE) = soil
  (; θ_sat) = soil.param
  zwt = soil.zwt * 100 # [m] -> [cm]
  z = soil.z_cm

  z[N+1] = 0.5 * (zwt + z[N])
  Δz[N+1] = zwt >= z[N] ? Δz[N] : abs(zwt - z[N])
  cal_θEψE!(soil) # update θE, ψE

  for i in 1:N
    i2 = min(N, i + 1)
    _θ = 0.5(θ[i] + θ[i2])
    _θsat = 0.5(θ_sat[i] + θ_sat[i2])
    se = clamp(_θ / _θsat, 0.01, 1.0)
    ψ[i] = Retention_ψ_Se(se, param[i])
  end

  i = N
  if jwt == N # GW under soil profile
    se = 0.5 * (θ_sat[i] + θ[i]) / θ_sat[i]
    se = clamp(se, 0.01, 1.0)
    ψ[i+1] = Retention_ψ_Se(se, param[i]) # N+1层的ψ，用的是第N层
  end
  jwt < N && (ψ[i+1] = 0.0)

  ## kernal
  for i = 1:N
    dz = (z[i+1] - z[i])
    dψ = (ψ[i+1] - ψE[i+1]) - (ψ[i] - ψE[i])
    Q[i] = -K₊ₕ[i] * dψ / dz
  end
  jwt < N && (Q[N] = 0.0)
  return Q0
end


function RichardsEquation_Zeng2009(dθ::AbstractVector{T}, θ::AbstractVector{T}, p::Soil{T}, t; method="ψ0") where {T<:Real}
  p.timestep += 1
  (; ibeg, N, Q, Q0, sink) = p # Δz, z, 
  Δz = p.Δz_cm
  Q0 = cal_Q_Zeng2009!(p, θ; Q0)

  dθ[ibeg] = ((-Q0 + Q[ibeg]) - sink[ibeg]) / Δz[ibeg] / 3600.0
  @inbounds for i in ibeg+1:N
    dθ[i] = ((-Q[i-1] + Q[i]) - sink[i]) / Δz[i] / 3600.0 # [m3 m-3] / h-1
  end
end

export cal_Q_Zeng2009!, RichardsEquation_Zeng2009
