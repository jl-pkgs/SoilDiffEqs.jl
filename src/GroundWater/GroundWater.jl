using UnPack
export find_jwt, GW_Rsb
export GW_Update_ZWT!, GW_Correctθ!

include("GW_Update_ZWT.jl")
include("GW_Correctθ.jl")
include("specific_yield.jl")

# 水位向下为正，地表为0
# 0 ~ N + 1
function find_jwt(z₊ₕ::AbstractVector, zwt::Real; N::Integer=length(z₊ₕ))
  # N = length(z₊ₕ)
  zwt > 0 && return 0
  zwt == 0 && return 1
  zwt < z₊ₕ[N] && return N + 1

  for j in 1:N
    zwt >= z₊ₕ[j] && return j
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

