# https://chatgpt.com/c/67aad86c-8458-8012-85e7-07d3ae77aa16
import HypergeometricFunctions: _₂F₁

# z: 向下为负
function _cal_θE_campbell(z1::T, z0::T, zwt::T, ψ_sat::T, par::ParamCampbell{T}) where {T<:Real}
  (; θ_sat, b) = par
  C = ψ_sat + zwt
  Δz = z0 - z1
  u1 = ((C - z1) / ψ_sat)^(1 - 1 / b) # z1 = z[j]
  u0 = ((C - z0) / ψ_sat)^(1 - 1 / b) # z0 = z[j-1]
  θE = ψ_sat * θ_sat / (1 - 1 / b) / Δz * (u1 - u0) # Zeng2009, Eq.9
  θE = clamp(θE, 0, θ_sat)
  return θE
end

function cal_θE(z1::T, z0::T, zwt::T, ψ_sat::T, par::ParamCampbell{T}) where {T<:Real}
  (; θ_sat) = par
  if zwt >= z0
    θE = θ_sat
  elseif z0 < zwt < z1
    d1 = zwt - z0 # 未饱和
    d2 = z1 - zwt   # 饱和
    _θE = _cal_θE_campbell(zwt, z0, zwt, ψ_sat, par)
    _θE = (_θE * d1 + θ_sat * d2) / (d1 + d2)
    θE = clamp(_θE, 0, θ_sat)
  elseif zwt <= z1
    θE = _cal_θE_campbell(z1, z0, zwt, ψ_sat, par)
  end
  return θE
end

## van Genuchten 1980
# 辅助函数 F(ψ) = ψ * hyp2f1(1/n, 1-1/n, 1+1/n, -(-αψ)^n)
function F_van1980(ψ::T, α::T, n::T) where {T<:Real}
  αψ = abs(-α * ψ)
  return ψ * _₂F₁(1 / n, 1 - 1 / n, 1 + 1 / n, -(αψ)^n)
end

# 解析积分函数，利用超几何函数计算积分的解析表达式
# 公式为：
#   I = θ_r*(z_{i+1/2}-z_{i-1/2})
#     + (θ_s-θ_r)*[F(C-z_{i-1/2}) - F(C-z_{i+1/2})]
function _cal_θE_van1980(z1::T, z0::T, zwt::T, ψ_sat::T, par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, α, n) = par
  Δz = z0 - z1 # make sure positive
  C = ψ_sat + zwt
  I1 = θ_res * (z1 - z0)
  I2 = (θ_sat - θ_res) * (F_van1980(C - z0, α, n) - F_van1980(C - z1, α, n))
  return -(I1 + I2) / Δz
end

function cal_θE(z1::T, z0::T, zwt::T, ψ_sat::T,
  par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat) = par
  
  if zwt >= z0
    θE = θ_sat
  elseif z1 <= zwt < z0
    d1 = zwt - z0 # 未饱和
    d2 = z1 - zwt   # 饱和
    _θE = _cal_θE_van1980(zwt, z0, zwt, ψ_sat, par)
    _θE = (_θE * d1 + θ_sat * d2) / (d1 + d2)
    θE = clamp(_θE, 0, θ_sat)
  elseif zwt < z1
    θE = _cal_θE_van1980(z1, z0, zwt, ψ_sat, par)
  end
  return θE
end

export cal_θE
