## 参考CoLM的方案，考虑地下水的补给与排泄，只考虑地下水垂直方向与土壤水的交互作用。
using SoilDifferentialEquations

# 从下至上，不受地下水影响的第一层
function find_jwt(soil::Soil{T}, zwt::T) where {T<:Real}
  (; N, z₊ₕ) = soil # 注意z是负值
  jwt = 0
  for i = 1:N
    if zwt >= z₊ₕ[i]
      jwt = i - 1
      break
    end
  end
  return jwt
end


using Test

@testset "find_jwt" begin
  @test find_jwt(soil, -0.01) == 0
  @test find_jwt(soil, -0.02) == 0
  @test find_jwt(soil, -0.03) == 1
end

N = 100
dz = fill(0.02, N)
θ = fill(0.3, N)
soil = Soil(dz; θ)

# 地下水的补给
(; θ, wa) = soil
zwt = -0.05
# Sy[1] = 2.0

# Qr: recharge rate，注意是正值
# Qd: discharge rate
function GW_Updateθ!(soil::Soil{T}, Δt, Qr, Qd) where {T<:Real}
  (; Sy, zwt, z₊ₕ, N) = soil
  Sy = specific_yield!(soil, zwt)

  ∑Qr = Qr * Δt
  jwt = find_jwt(soil, zwt) # 不受地下水影响的第一层

  if ∑Qr > 0.0
    # 水位上升
    for i = jwt+1:N
      _sy = Sy[i] # unitless
      _recharge = clamp(∑Qr, 0, -_sy * (zwt - z₊ₕ[i-1]) * Δt) # 补给量，正值
      zwt = min(zwt + _recharge / _sy, 0) # 如果补给增加，则水位应该上升，水位最大为0
      ∑Qr -= _recharge
      ∑Qr <= 0 && break
    end
  else
    # 水位下降
    for i = jwt+1:N
      _sy = Sy[i] # unitless
      _discharge = clamp(∑Qr, _sy * (z₊ₕ[i] - zwt) * Δt, 0) # 排泄量，负值

      ∑Qr -= _discharge
      if ∑Qr > 0
        zwt = min(zwt + _discharge / _sy, 0) # 加一个负值，水位下降
        break
      else
        zwt = z₊ₕ[i] # 否则这一层排空
      end
    end
    # 放到最后一层，不可能出现这种情况
    ∑Qr > 0 && (zwt = zwt + ∑Qr/Sy[N]) 
  end
end

