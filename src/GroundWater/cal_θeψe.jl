# Constants
const SHR_CONST_TKFRZ = 273.15
const SHR_CONST_LATICE = 3.337e5
const SHR_CONST_G = 9.80665


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
function cal_θeψe!(θE, ψE, z, zwt, jwt; θ_sat, ψ_sat, B, ψmin)
  N = length(z)
  for j = 1:N
    z0 = j == 1 ? 0 : z[j-1]
    z1 = z[j]
    θE[j]= cal_θE(z1, z0, zwt, ψ_sat[j], θ_sat[j], B[j]; use_clamp)
    ψE[j] = cal_ψ(θE[j], θ_sat[j], ψ_sat[j], B[j]; ψmin)
  end

  # If water table is below soil column calculate ψE for the 11th layer
  if jwt == N
    # 积分的过程在：z₊ₕ[N] ~ zwt，因此，需要注意，最后一层的θ_E、ψ_E代表的区间
    j = N
    θE[j+1] = cal_θE(zwt, z[j], ψ_sat[j], zwt, B[j])
    ψE[j+1] = cal_ψ(θE[j+1], θ_sat[j], ψ_sat[j], B[j]; ψmin)
  end
  ψE
end
