"""
- `Δt`       : [h]
- `drainage` : [cm h-1], 排泄为正
- `wa`       : [mm]
- `zwt`      : [m]
"""
function GW_Correctθ!(soil::Soil{FT,P}, θ::AbstractVector{FT}, zwt, wa, Δt, drainage) where {FT<:Real,P}
  (; N, Δz, Sy_e) = soil
  (; θ_sat) = soil.param

  # jwt = find_jwt(z₊ₕ, zwt)
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
        uex = exceed * 1000 # [m] to [cm]
      else
        θ[j-1] = θ[j-1] + exceed / Δz[j-1]
      end
    end
  end

  ## 2. 亏损？
  ∑_neg = 0.0
  for j = 1:N
    if θ[j] < 0
      ∑_neg += θ[j] * Δz[j] # [m]
      θ[j] = 0.0
    end
  end

  # 1. 少排一点水; 2. drainage扣完之后，从地下水中扣除
  drainage += ∑_neg / Δt * 100 # [m h-1] to [cm h-1]
  # @show ∑_neg, drainage
  if drainage < 0
    wa += drainage * Δt * 10 # [cm] to [mm]
    zwt += drainage * Δt / Sy_e[N] / 100
    drainage = 0.0
  end
  (; wa, uex, drainage)
end
