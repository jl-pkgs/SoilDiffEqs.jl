"""
> 按照CoLM的方式进行修改

Kun的做法是合理的，非饱和带和饱和带应该分开，否则地下水层会出现镂空。

该函数仅处理非饱和带的土壤移动，非饱和带土壤水分的变化不会引起地下水水位的变化。
饱和带的水分移动，会引起地下水水位的变化。

# TODO: 
1. 补充一个均衡水位的示意图
"""

"""
    GW_Update_ZWT!(soil::Soil, θ::AbstractVector, zwt, wa, ∑;)

> This function is for `SoilDiffEqs`, kdd, 20250516
目的在于处理QN和地下径流引起的地下水水位的上升或下降。

- `∑`: `∑ = drainage * dt / 100`, [cm h-1 * h-1 / 100] = [cm/100] = [m]
  > `∑ = -Q[N] * 1h / 100`, in [m]
  + `∑ > 0`，SM补给GW，GW上升
  + `∑ < 0`，GW补给SM，GW下降
"""
function GW_Update_ZWT!(soil::Soil, θ::AbstractVector, zwt, wa, ∑;) #where {T<:Real}
  (; z₊ₕ, Δz, N) = soil
  (; θ_sat, θ_fc) = soil.param

  j = find_jwt(z₊ₕ, zwt)
  # _θ_wp = 0.01
  # _θ_wp = i == 1 ? 0.01 : θ_wp # 土壤蒸发可超越凋萎含水量的限制

  ## 非饱和带
  if ∑ >= 0 # 补给、水位上升
    for i = j:-1:1
      _sy = i == N + 1 ? Sy[N] : Sy[i] # unitless
      z0 = i == 1 ? 0.0 : z₊ₕ[i-1]
      z1 = i == N+1 ? zwt : z₊ₕ[i]

      _recharge = clamp(∑, 0, _sy * (z0 - z1)) # 补给量，正值
      zwt = zwt + _recharge / _sy 
      # 如果补给增加，则水位应该上升，无需限制zwt，因为上界是_z₊ₕ

      i <= N && (θ[i] += _recharge / Δz[i]) # [cm] to [m^3/m^3]
      ∑ -= _recharge
      ∑ <= 0 && break
    end
    ∑ > 0 && (uex = ∑*1000) # excess water to soil surface
  end

  if ∑ < 0 # 排泄、水位下降
    for i = j:-1:1
      _sy = i == N + 1 ? Sy[N] : Sy[i] # unitless
      z0 = i == 1 ? 0.0 : z₊ₕ[i-1]  
      z1 = i == N + 1 ? zwt : z₊ₕ[i]

      _drainage = clamp(∑, -(θ_sat[i] - _sy) * (z0 - z1), 0) # 排泄量，负值
      zwt += _drainage / _sy

      i <= N && (θ[i] += _drainage / Δz[i]) # [cm] to [m^3/m^3]
      ∑ -= _drainage
      ∑ >= 0 && break
    end
    ∑ < 0 && (zwt += ∑ / (θ_sat[N] - Sy[N])) # excess water to soil surface
  end
  wa += ∑*1000 # in cm
  (; zwt, wa, uex)
end
