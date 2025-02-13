"""
  计算每一层的给水度

Sy = θ_sat - θ(ψ_sat + zwt)

前提条件是：每层土壤均接近饱和，才能得到这个公式，`ψ + z = ψ_sat + zwt`，即Q = 0。
"""
function specific_yield!(soil::Soil{T}, zwt::T; sy_max::T=0.02) where {T<:Real}
  (; N, Sy) = soil
  (; θ_sat, θ_res, α, n, m, use_m, b, method_retention) = soil.param

  if method_retention == "van_Genuchten"
    _ψ_sat = 0.0 # !Note
    for i = 1:N
      _m = use_m ? m[i] : 1 - 1 / n[i] # 不参与参数优化
      Sy[i] = θ_sat[i] - van_Genuchten_θ(_ψ_sat + zwt, θ_sat[i], θ_res[i], α[i], n[i], _m)
    end
  elseif method_retention == "Campbell"
    for i = 1:N
      Sy[i] = θ_sat[i] - Campbell_θ(ψ_sat[i] + zwt, ψ_sat[i], θ_sat[i], b[i])
    end
  end
  clamp!(Sy, 0, sy_max)
end


function cal_θKCap!(soil::Soil{T}, ψ::AbstractVector{T}) where {T<:Real}
  (; ibeg, N, θ, K, Cap) = soil
  (; θ_sat, θ_res, ψ_sat, Ksat, α, n, m, use_m, b, method_retention) = soil.param

  if method_retention == "van_Genuchten"
    @inbounds for i in ibeg:N
      # _m = m[i]
      _m = use_m ? m[i] : 1 - 1 / n[i] # 不参与参数优化
      θ[i], K[i], Cap[i] = van_Genuchten(ψ[i], θ_sat[i], θ_res[i], Ksat[i], α[i], n[i], _m)
    end
  elseif method_retention == "Campbell"
    @inbounds for i in ibeg:N
      θ[i], K[i], Cap[i] = Campbell(ψ[i], ψ_sat[i], θ_sat[i], Ksat[i], b[i])
    end
  end
end

function K_interface(K1::T, K2::T, z1::T, z2::T) where {T<:Real}
  # K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  K1 * K2 * (z1 + z2) / (K1 * z2 + K2 * z1) # Eq. 5.16, 
end

function cal_K₊ₕ!(soil::Soil)
  (; N, ibeg, Δz, K, K₊ₕ) = soil
  
  @inbounds for i = ibeg:N-1
    # d1 = z[i] - z₊ₕ[i]
    # d2 = z₊ₕ[i] - z[i+1]
    # K₊ₕ[i] = K[i] * K[i+1] * (d1 + d2) / (K[i] * d2 + K[i+1] * d1) 
    K₊ₕ[i] = K_interface(K[i], K[i+1], Δz[i], Δz[i+1])
  end
end

cal_K!(soil::Soil) = cal_K!(soil, soil.θ)
function cal_K!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, K) = soil
  (; θ_sat, θ_res, Ksat, m, b, method_retention) = soil.param

  if method_retention == "van_Genuchten"
    @inbounds for i = ibeg:N
      K[i] = van_Genuchten_K(θ[i], θ_sat[i], θ_res[i], Ksat[i], m[i])
    end
  elseif method_retention == "Campbell"
    @inbounds for i = ibeg:N
      K[i] = Campbell_K(θ[i], θ_sat[i], Ksat[i], b[i])
    end
  end
end


cal_ψ!(soil::Soil) = cal_ψ!(soil, soil.θ)
function cal_ψ!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, ψ) = soil
  (; θ_sat, θ_res, α, n, m, ψ_sat, b, method_retention) = soil.param
  if method_retention == "van_Genuchten"
    @inbounds for i = ibeg:N
      ψ[i] = van_Genuchten_ψ(θ[i], θ_sat[i], θ_res[i], α[i], n[i], m[i])
    end
  elseif method_retention == "Campbell"
    @inbounds for i = ibeg:N
      ψ[i] = Campbell_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
    end
  end
end


function Init_ψ0(soil::Soil{T}, θ::T) where {T<:Real}
  i = soil.ibeg # 默认采用第一层进行初始化
  (; θ_sat, θ_res, α, n, m, ψ_sat, b, method_retention) = soil.param

  if method_retention ==
    "van_Genuchten"
    ψ0 = van_Genuchten_ψ(θ, θ_sat[i], θ_res[i], α[i], n[i], m[i])
  elseif method_retention == "Campbell"
    ψ0 = Campbell_ψ(θ, θ_sat[i], ψ_sat[i], b[i])
  end
  return ψ0
end


export specific_yield!
export cal_K!, cal_ψ!, cal_θKCap!, cal_K₊ₕ!
export Init_SoilWaterParam, Init_ψ0
