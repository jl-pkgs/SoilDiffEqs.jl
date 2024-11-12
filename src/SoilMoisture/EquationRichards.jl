function soil_WaterFlux!(soil::Soil{T}, θ::AbstractVector{T};
  ψ0::T=NaN, Q0::T=NaN, method="ψ0") where {T<:Real}

  (; ibeg, N, Q, K, K₊ₕ, ψ) = soil
  z = soil.z_cm

  # need to update here
  cal_K!(soil, θ)
  cal_ψ!(soil, θ)

  if method == "ψ0"
    z_prev = ibeg == 1 ? 0 : z[ibeg-1]
    Q0 = -K[ibeg] * ((ψ0 - ψ[ibeg]) / (z_prev - z[ibeg]) + 1) # [cm/s]
  elseif method == "Q0"
    # Q0 = Q0
  end

  cal_K₊ₕ!(soil)
  @inbounds for i in ibeg:N-1
    # K₊ₕ = (K[i] + K[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    Q[i] = -K₊ₕ[i] * ((ψ[i] - ψ[i+1]) / Δz₊ₕ + 1)
  end
  Q[N] = -K[N] # 尾部重力排水
  Q0
end

"Update soil water content θ according to Q"
function soil_Updateθ!(soil::Soil{T}, dt::Real; method="ψ0") where {T<:Real}
  (; ibeg, N, θ, Q, ψ0, Q0, sink) = soil # Δz, z, 
  Δz = soil.Δz_cm                               # [cm]
  Q0 = soil_WaterFlux!(soil, θ; ψ0, Q0, method) # [cm/s]
  (; θ_sat, θ_res) = soil.param

  # TODO: 
  # 若某一层发生了饱和，则继续向下传导
  # 若某一层的土壤水分全部耗干，则继续向下抽水
  θ[ibeg] += -((Q0 - Q[ibeg]) + sink[ibeg]) * dt / Δz[ibeg]
  θ[ibeg] = clamp(θ[ibeg], θ_res[ibeg], θ_sat[ibeg])

  @inbounds for i in ibeg+1:N
    θ[i] += -((Q[i-1] - Q[i]) + sink[i]) * dt / Δz[i] # [m3 m-3]
    θ[i] = clamp(θ[i], θ_res[i], θ_sat[i])
  end
end

"""
# Arguments
- `method`: 
  + `ψ0`: ψ0 boundary condition, 第一类边界条件
  + `Q0`: Q0 boundary condition, 第二类边界条件
"""
function RichardsEquation(dθ::AbstractVector{T}, θ::AbstractVector{T}, p::Soil{T}, t; method="ψ0") where {T<:Real}
  p.timestep += 1
  (; ibeg, N, Q, ψ0, Q0, sink) = p # Δz, z, 
  Δz = p.Δz_cm
  Q0 = soil_WaterFlux!(p, θ; ψ0, Q0, method)

  dθ[ibeg] = -((Q0 - Q[ibeg]) + sink[ibeg]) / Δz[ibeg]
  @inbounds for i in ibeg+1:N
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
