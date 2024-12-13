# Constants
const SHR_CONST_TKFRZ = 273.15
const SHR_CONST_LATICE = 3.337e5
const SHR_CONST_G = 9.80665


function cal_θE(z1::T, z0::T, zwt::T, ψ_sat::T, B::T; use_clamp=true) where {T<:Real}
  C = ψ_sat + zwt
  Δz = z1 - z0
  u1 = ((C - z1) / ψ_sat)^(1 - 1 / B) # z1 = z[j]
  u0 = ((C - z0) / ψ_sat)^(1 - 1 / B) # z2 = z[j-1]
  θE = -ψ_sat[j] * θ_sat[j] / (1 - 1 / B) / Δz * (u1 - u0) # Zeng2009, Eq.9
  use_clamp && (θE = clamp(θE, 0, θ_sat))
  return θE
end

function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, B::T; ψmin::T) where {T<:Real}
  se = max(θ / θ_sat, 0.01, 1.0)
  ψ = -ψ_sat * se^(-B)
  return max(ψ, ψmin)
end

function cal_ψ(se::T, ψ_sat::T, B::T; ψmin::T) where {T<:Real}
  # se = max(θ / θ_sat, 0.01, 1.0)
  ψ = -ψ_sat * se^(-B)
  return max(ψ, ψmin)
end


# Calculate the equilibrium water content based on the water table depth
function cal_θeψe!(θE, ψE, z, zwt, jwt; θ_sat, ψ_sat, bsw, ψmin)
  N = length(z)
  for j = 1:N
    z0 = j == 1 ? 0 : z[j-1]
    z1 = z[j]

    if zwt <= z0
      # 全部饱和
      θE[j] = θ_sat[j]
    elseif z0 < zwt < z1
      # 部分饱和
      d1 = zwt - z0 # 未饱和
      d2 = z1 - zwt   # 饱和
      _θE = cal_θE(zwt, z0, zwt, ψ_sat[j], bsw[j]; use_clamp=false)
      θE[j] = (_θE * d1 + θ_sat[j] * d2) / (d1 + d2)
      θE[j] = clamp(θE[j], 0, θ_sat[j])
    else
      # 非饱和
      θE[j] = cal_θE(z1, z0, zwt, ψ_sat[j], bsw[j])
    end
    ψE[j] = cal_ψ(θE[j], θ_sat[j], ψ_sat[j], bsw[j]; ψmin)
  end

  # If water table is below soil column calculate ψE for the 11th layer
  if jwt == N
    # 积分的过程在：z₊ₕ[N] ~ zwt，因此，需要注意，最后一层的θ_E、ψ_E代表的区间
    j = N
    θE[j+1] = cal_θE(zwt, z[j], ψ_sat[j], zwt, bsw[j])
    ψE[j+1] = cal_ψ(θE[j+1], θ_sat[j], ψ_sat[j], bsw[j]; ψmin)
  end
  ψE
end
