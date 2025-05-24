"""
- Q0: Infiltration rate, [cm h-1], should be negative or zero
"""
function cal_Q_Zeng2009!(soil::Soil{T}, θ::AbstractVector{T}) where {T<:Real}
  (; N, jwt, Q, ψ, ψE, K₊ₕ) = soil
  (; θ_sat, θ_res, param) = soil.param
  zwt = soil.zwt * 100 # [m] -> [cm]
  z = soil.z_cm
  Δz = soil.Δz_cm

  cal_θEψE!(soil) # update θE, ψE
  cal_K!(soil, θ)
  cal_ψ!(soil, θ)

  z[N+1] = 0.5 * (zwt + z[N])
  Δz[N+1] = zwt >= z[N] ? Δz[N] : abs(zwt - z[N])

  i = N
  if jwt == N # GW under soil profile
    _θ = 0.5 * (θ[i] + θ_sat[i])
    se = clamp((_θ - θ_res[i]) / (θ_sat[i] - θ_res[i]), 0.01, 1.0)
    # se = 0.5 * (θ_sat[i] + θ[i]) / θ_sat[i]
    # se = clamp(se, 0.01, 1.0)
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
  return Q
end


function RichardsEquation_Zeng2009(dθ::AbstractVector{T}, θ::AbstractVector{T}, soil::Soil{T}, t; method="ψ0") where {T<:Real}
  # 这里采用的是Q0
  soil.timestep += 1
  (; ibeg, N, Q, Q0, sink) = soil # Δz, z, 
  Δz = soil.Δz_cm
  cal_Q_Zeng2009!(soil, θ)

  dθ[ibeg] = ((-Q0 + Q[ibeg]) - sink[ibeg]) / Δz[ibeg] / 3600.0
  @inbounds for i in ibeg+1:N
    dθ[i] = ((-Q[i-1] + Q[i]) - sink[i]) / Δz[i] / 3600.0 # [m3 m-3] / h-1
  end
end


export cal_Q_Zeng2009!, RichardsEquation_Zeng2009

