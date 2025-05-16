"""
    GW_Update_ZWT!(soil, soilpar; Q_in=0.0)

> This function is for `SoilDiffEqs`, kdd, 20250516
目的在于处理QN和地下径流引起的地下水水位的上升或下降。

> 按照CoLM的方式进行修改

Kun的做法是合理的，非饱和带和饱和带应该分开，否则地下水层会出现镂空。

该函数仅处理非饱和带的土壤移动，非饱和带土壤水分的变化不会引起地下水水位的变化。
饱和带的水分移动，会引起地下水水位的变化。

# TODO: 
1. 补充一个均衡水位的示意图
"""
function GW_Update_ZWT!(soil::Soil; ) #where {T<:Real}
  (; z₊ₕ, Δz, zwt, θ, Q, N) = soil
  (; θ_sat, θ_wp) = soil.param

  # 注意单位
  j = find_jwt(z₊ₕ, zwt)
  ∑ = -Q[N]
  # `∑ > 0`，SM补给GW，GW上升
  # `∑ < 0`，GW补给SM，GW下降
  # _θ_wp = 0.01
  # _θ_wp = i == 1 ? 0.01 : θ_wp # 土壤蒸发可超越凋萎含水量的限制

  ## 非饱和带
  if ∑ >= 0
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
    ∑ > 0 && (uex = ∑) # excess water to soil surface
  end

  ## 补给和排泄，都是从最下层开始
  if ∑ < 0
    for i = j:N
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
  wa += ∑ # in cm
  (; zwt, wa, uex)
end



# drainage < 0, [m s-1]
function GW_Correctθ!(soil::Soil{FT,P}, θ::AbstractVector{FT}, zwt, wa, Δt, drainage) where {FT<:Real,P}
  (; N, Δz, z₊ₕ, Sy) = soil
  (; θ_sat) = soil.param

  jwt = find_jwt(z₊ₕ, zwt)
  zwt = clamp(zwt, 0.0, 80.0) # 地下水水位在[0, 80m]
  # 强制限制水位，不考虑水量平衡是否合适？

  uex = 0.0
  ## 1. 超饱和？
  exceed = 0.0
  for j = N:-1:1
    exceed = max((θ[j] - θ_sat[j]) * Δz[j], 0.0)
    if exceed > 0.0
      θ[j] = θ_sat[j]
      if j == 1
        uex = exceed
      else
        θ[j-1] = θ[j-1] + exceed / Δz[j-1]
      end
    end
  end

  ## 2. 亏损？
  ∑_neg = 0.0
  for j = 1:N
    if θ[j] < 0
      ∑_neg += θ[j] * Δz[j]
      θ[j] = 0.0
    end
  end

  # 1. 少排一点水; 2. drainage扣完之后，从地下水中扣除
  drainage = drainage + ∑_neg / Δt * 1000 # [mm s-1] or [kg s-1]
  if drainage < 0
    wa += drainage * Δt
    zwt = zwt + drainage * Δt / Sy[N] / 1000
    drainage = 0.0
  end
  (; wa, uex, drainage)
end
