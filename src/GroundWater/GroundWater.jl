using UnPack
export find_jwt, GW_Rsb
export GW_Update_ZWT!, GW_Correctθ!


"""
  计算每一层的给水度

Sy = θ_sat - θ(ψ_sat + zwt) # 可以理解为孔隙度

前提条件是：每层土壤均接近饱和，才能得到这个公式，`ψ + z = ψ_sat + zwt`，即Q = 0。
"""
function specific_yield!(soil::Soil{T}, zwt::T; sy_max::T=0.02) where {T<:Real}
  (; N, Sy) = soil
  (; θ_sat, param, method_retention) = soil.param

  iszero_ψ = method_retention == "van_Genuchten"
  for i = 1:N
    _ψ_sat = iszero_ψ ? 0.0 : θ_sat[i] # !Note
    Sy[i] = θ_sat[i] - Retention_θ(_ψ_sat + zwt, param[i]) # 这里存在更好的计算方法，采用θE
  end
  clamp!(Sy, 0.0, sy_max)
end

# 水位向下为正，地表为0
# 0 ~ N + 1
function find_jwt(z₊ₕ::AbstractVector, zwt::Real)
  N = length(z₊ₕ)
  zwt <= 0 && return 0
  zwt >= z₊ₕ[end] && return N + 1

  for j in 1:N
    zwt <= z₊ₕ[j] && return j
  end
end

# function find_jwt(z₊ₕ::AbstractVector, zwt::Real; N::Int=length(z₊ₕ))
#   jwt = N
#   zwt_abs = abs(zwt)
#   for j in 1:N
#     if zwt_abs <= abs(z₊ₕ[j])
#       jwt = j - 1
#       break
#     end
#   end
#   return jwt
# end

function GW_Rsb(zwt::Real)
  R_sb_max = 39.0 # mm day-1
  f = 1.25e-3     # mm-1
  return R_sb_max * exp(-f * zwt)
end


include("GW_Update_ZWT.jl")
