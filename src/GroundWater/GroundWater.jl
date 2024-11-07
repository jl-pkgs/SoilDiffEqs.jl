using UnPack
export find_jwt, GW_UpdateRecharge!, GW_UpdateDrainage!, GW_Correctθ!

# 水位为正
function find_jwt(z₊ₕ::AbstractVector, zwt::Real)
  N = length(z₊ₕ)
  jwt = N

  zwt_abs = abs(zwt)
  for j in 1:N
    if zwt_abs <= abs(z₊ₕ[j])
      jwt = j - 1
      break
    end
  end
  return jwt
end


include("GW_Update.jl")
