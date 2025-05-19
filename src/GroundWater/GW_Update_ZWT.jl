"""
    GW_Update_ZWT!(soil::Soil, θ::AbstractVector, zwt, wa, ∑;)

> This function is for `SoilDiffEqs`, kdd, 20250516
目的在于处理QN和地下径流引起的地下水水位的上升或下降。

## Arguments
- `wa`: wa定义为地下水水量，不包含土壤水的部分

- `∑` : [mm], `∑ = drainage * dt * 10`, [cm h-1 * h *10] = [mm]
  > `∑ = -Q[N] * 1h * 10`, in [mm]
  + `∑ > 0`，SM补给GW，GW上升
  + `∑ < 0`，GW补给SM，GW下降
"""
function GW_Update_ZWT!(soil::Soil, θ::AbstractVector, zwt, wa, ∑;) #where {T<:Real}
  specific_yield!(soil, zwt)
  (; z₊ₕ, Δz, N, Sy_r, Sy_d, Sy_e) = soil
  (; θ_sat) = soil.param
  wa += ∑ # [mm]
  Sy = ∑ >= 0 ? Sy_r : Sy_d
  j = find_jwt(z₊ₕ, zwt; N)
  # _θ_wp = 0.01
  # _θ_wp = i == 1 ? 0.01 : θ_wp # 土壤蒸发可超越凋萎含水量的限制

  uex = 0.0
  ## 非饱和带
  if ∑ >= 0 # 补给、水位上升
    for i = j:-1:1
      # _sy = Sy[i] # unitless
      _i = min(i, N)
      _θsat = θ_sat[_i]
      z0 = i == 1 ? 0.0 : z₊ₕ[i-1]
      z1 = i == N + 1 ? zwt : z₊ₕ[i]
      d = z0 - z1

      _θ = i == N + 1 ? (θ[N] + _θsat) / 2 : θ[i]
      d1 = z0 - zwt # 未饱和
      d2 = zwt - z1 # 饱和，饱和部分排泄
      # θ_unsat * d1 + θ_sat * d2 = θ[i] * d
      θ_unsat = (_θ * d - _θsat * d2) / d1 # 如果为负，要如何处理？
      _sy = _θsat - θ_unsat

      _recharge = clamp(∑ / 1000, 0, _sy * d1) # 补给量，正值
      zwt = zwt + _recharge / _sy
      # 无需限制zwt，因为上界是_z₊ₕ

      i <= N && (θ[i] += _recharge / Δz[i]) # [cm] to [m^3/m^3]
      ∑ -= _recharge * 1000
      ∑ <= 0 && break
    end
    ∑ > 0 && (uex = ∑) # excess water to soil surface
  end

  if ∑ < 0 # 排泄、水位下降
    for i = j:N
      _sy = Sy[i] # unitless
      z0 = i == 1 ? 0.0 : z₊ₕ[i] # zwt
      z1 = z₊ₕ[i]
      # d1 = z0 - zwt # 未饱和
      d2 = zwt - z1 # 饱和，饱和部分排泄

      _drainage = clamp(∑ / 1000, -_sy * d2, 0) # 排泄量，负值
      zwt += _drainage / _sy # 向下降d2, 到z₊ₕ[i]

      θ[i] += _drainage / Δz[i] # [cm] to [m^3/m^3]
      ∑ -= _drainage * 1000
      ∑ >= 0 && break
    end
    ∑ < 0 && (zwt += (∑ / 1000) / Sy[N]) # excess water to soil surface
  end
  (; zwt, wa, uex)
end


GW_Update_ZWT!(soil::Soil, θ::AbstractVector; zwt, wa, ∑) =
  GW_Update_ZWT!(soil, θ, zwt, wa, ∑)
