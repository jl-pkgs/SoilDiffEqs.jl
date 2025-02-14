export cal_Q!

function cal_Q!(soil::Soil{T}; Q0::T=0.0, ψ0::T=NaN)::T where {T<:Real}
  (; θ) = soil
  method = "Q0"
  !isnan(ψ0) && (method = "ψ0")
  Q0 = cal_Q!(soil, θ; Q0, ψ0, method)
  return Q0
end

function cal_Q!(soil::Soil{T}, θ::AbstractVector{T};
  ψ0::T=NaN, Q0::T=NaN, method="ψ0")::T where {T<:Real}

  (; ibeg, N, Q, K, K₊ₕ, ψ) = soil
  z = soil.z_cm
  Δz₊ₕ = soil.Δz₊ₕ_cm
  # need to update here
  cal_K!(soil, θ)
  cal_ψ!(soil, θ)

  if method == "ψ0"
    # _dz::T = z[ibeg-1] - z[ibeg]
    _dz::T = ibeg == 1 ? -z[1] : Δz₊ₕ[ibeg-1]
    _K₊ₕ::T = ibeg == 1 ? K[1] : K₊ₕ[ibeg-1]
    Q0 = -_K₊ₕ * ((ψ0 - ψ[ibeg]) / _dz + 1) # [cm h-1]
  end

  @inbounds for i in ibeg:N-1
    Q[i] = -K₊ₕ[i] * ((ψ[i] - ψ[i+1]) / Δz₊ₕ[i] + 1.0)
  end
  Q[N] = -K[N] # 尾部重力排水
  soil.Q0 = Q0
  return Q0
end

# dt: in hours
"Update soil water content θ according to Q"
function soil_Updateθ!(soil::Soil{T}, dt::T) where {T<:Real}
  (; ibeg, N, θ, Q, Q0, sink) = soil # Δz, z, 
  (; θ_sat, θ_res) = soil.param
  Δz = soil.Δz_cm                               # [cm]
  dt = dt / 3600.0 # [s] -> [h]
  # Q0 = cal_Q!(soil, θ; ψ0, Q0, method) # [cm h-1]

  # TODO: 
  # 若某一层发生了饱和，则继续向下传导
  # 若某一层的土壤水分全部耗干，则继续向下抽水
  if ibeg > 1
    θ[ibeg-1] -= -Q0 * dt / Δz[ibeg-1]
  end
  θ[ibeg] += ((-Q0 + Q[ibeg]) - sink[ibeg]) * dt / Δz[ibeg]
  # θ[ibeg] = clamp(θ[ibeg], θ_res[ibeg], θ_sat[ibeg])

  @inbounds for i in ibeg+1:N
    θ[i] += ((-Q[i-1] + Q[i]) - sink[i]) * dt / Δz[i] # [m3 m-3]
    # θ[i] = clamp(θ[i], θ_res[i], θ_sat[i])
  end
end

"""
# Arguments
- `method`: 
  + `ψ0`: ψ0 boundary condition, 第一类边界条件
  + `Q0`: Q0 boundary condition, 第二类边界条件
> dθ: actually is dθ/dt
"""
function RichardsEquation(dθ::AbstractVector{T}, θ::AbstractVector{T}, p::Soil{T}, t; method="ψ0") where {T<:Real}
  p.timestep += 1
  # mod(p.timestep, 1000) == 0 && println("timestep = ", p.timestep)

  (; ibeg, N, Q, ψ0, Q0, sink) = p # Δz, z, 
  Δz = p.Δz_cm
  Q0 = cal_Q!(p, θ; ψ0, Q0, method)

  dθ[ibeg] = -((Q0 - Q[ibeg]) + sink[ibeg]) / Δz[ibeg] / 3600.0
  @inbounds for i in ibeg+1:N
    _dθ = -(Q[i-1] - Q[i]) / Δz[i] - sink[i] / Δz[i] # [m3 m-3] / h-1
    dθ[i] = _dθ / 3600.0 # [m3 m-3] / s-1
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
