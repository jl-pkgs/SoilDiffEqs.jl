using UnPack
export find_jwt, GW_Rsb
export GW_Update_ZWT!, GW_Correctθ!

include("GW_Update_ZWT.jl")
include("GW_Correctθ.jl")


"""
  计算每一层的给水度

Sy = θ_sat - θ(ψ_sat + zwt) # 可以理解为孔隙度

前提条件是：每层土壤均接近饱和，才能得到这个公式，`ψ + z = ψ_sat + zwt`，即Q = 0。
"""
function specific_yield!(soil::Soil{T}, zwt::T;) where {T<:Real}
  (; N, Sy) = soil
  (; θ_sat, param, method_retention) = soil.param
  cal_θEψE!(soil) # update θE, ψE
  
  for i = 1:N+1
    _i = min(i, N) # N+1层，采用第N层的参数
    _ψ_sat = ψ_sat[_i]

    # 排泄与补给不同
    Sy_d[i] = θ_sat[_i] - θ_fc[_i] # 排泄
    Sy_r[i] = θ_sat[_i] - θ[_i]    # 补给
    Sy_e[i] = θ_sat[_i] - Retention_θ(_ψ_sat + zwt, param[_i]) # 均衡状态
    # 详见: CoLM TechNote, 2024, ch12.1.2, P257
    # ψE + z = ψ_sat + zwt, θE = cal_ψ(ψE)
    # ψE = ψ_sat + zwt, (z = 0)
  end
end

# 水位向下为正，地表为0
# 0 ~ N + 1
function find_jwt(z₊ₕ::AbstractVector, zwt::Real; N::Integer=length(z₊ₕ))
  # N = length(z₊ₕ)
  zwt > 0 && return 0
  zwt < z₊ₕ[N] && return N + 1

  for j in 1:N
    zwt > z₊ₕ[j] && return j
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


