"""
- `Δt`       : [h]
- `drainage` : [cm h-1], 排泄为正
- `wa`       : [mm], 地下蓄水量
- `zwt`      : [m]
"""
function GW_Correctθ!(soil::Soil{FT,P}, θ::AbstractVector{FT}; zwt, exceed2surf=true) where {FT<:Real,P}
  (; N, Δz) = soil
  (; θ_sat) = soil.param

  jwt = find_jwt(soil.z₊ₕ, zwt; N)
  # zwt = clamp(zwt, -80.0, 0.0) # 地下水水位在[0, 80m]
  uex = 0.0
  exceed = 0.0
  ∑_exceed = 0.0 # [mm]
  ∑_deficit = 0.0 
  ## 1. 超饱和，则排出去；一般采用顶层排出，如果是底部排出则补给wa
  if exceed2surf
    # 1.1 超饱和的部分，作为地表径流
    for j = N:-1:1
      exceed = max((θ[j] - θ_sat[j]) * Δz[j], 0.0)
      if exceed > 0.0
        θ[j] = θ_sat[j]
        if j == 1
          uex = exceed * 1000 # [m] to [cm]
        else
          θ[j-1] += exceed / Δz[j-1]
        end
      end
    end
  else
    # 1.2 超饱和的部分，补给地下水
    for j = 1:N
      exceed = max((θ[j] - θ_sat[j]) * Δz[j], 0.0)
      if exceed > 0.0
        θ[j] = θ_sat[j]
        if j == N
          ∑_exceed = exceed * 1000 # [m] to [cm]
        else
          θ[j+1] += exceed / Δz[j+1]
        end
      end
    end
  end

  ## 2.1: fix `θ[j] < 0`
  for j = 1:N
    if θ[j] < 0
      ∑_deficit += θ[j] * Δz[j] * 1000 # [m]
      θ[j] = 0.0
    end
  end

  ## 2.2: 地下水以下应为饱和
  for j = jwt+1:N
    ## 判断饱和的部分和未饱和的部分
    if θ[j] < θ_sat[j]
      ∑_deficit -= (θ_sat[j] - θ[j]) * Δz[j] * 1000 # [m]
      θ[j] = θ_sat[j]
    end
  end

  ∑ = ∑_exceed + ∑_deficit # 底部的QN做一个微略的修正
  (; uex, ∑)
end


"""
1. Richards Equation
2. GW_Correctθ!
3. GW_Update_ZWT!
"""
