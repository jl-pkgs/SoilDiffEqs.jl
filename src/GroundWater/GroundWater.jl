# Constants
const SHR_CONST_TKFRZ = 273.15
const SHR_CONST_LATICE = 3.337e5
const SHR_CONST_G = 9.80665


# 水位为正
function find_jwt(z₊ₕ::AbstractVector, zwt::Real)
  N = length(z₊ₕ)
  jwt = N

  for j in 1:N
    if zwt <= z₊ₕ[j]
      jwt = j - 1
      break
    end
  end
  return jwt
end


function cal_θE(z1::T, z0::T, zwt::T, ψ_sat::T, B::T; use_clamp=true) where {T<:Real}
  C = ψ_sat + zwt
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
