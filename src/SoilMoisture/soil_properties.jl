"""
  计算每一层的给水度

Sy = θ_sat - θ(ψ_sat + zwt)

前提条件是：每层土壤均接近饱和，才能得到这个公式，`ψ + z = ψ_sat + zwt`，即Q = 0。
"""
function specific_yield!(soil::Soil{T}, zwt::T; sy_max::T=0.02) where {T<:Real}
  (; N, Sy) = soil
  (; θ_sat, ψ_sat, method_retention) = soil.param
  param = soil.param.param

  iszero_ψ = method_retention == "van_Genuchten"
  for i = 1:N
    _ψ_sat = iszero_ψ ? T(0) : ψ_sat[i]
    Sy[i] = θ_sat[i] - Retention_θ(_ψ_sat + zwt, param[i])
  end
  clamp!(Sy, 0, sy_max)
end


function cal_θKCap!(soil::Soil{T,P}, param::StructVector{P}, 
  ψ::AbstractVector{T}) where {T<:Real,P<:AbstractSoilParam{T}}
  (; ibeg, N, θ, K, Cap) = soil
  
  @inbounds for i in ibeg:N
    par = param[i]
    θ[i] = Retention_θ(ψ[i], par)
    K[i] = Retention_K(θ[i], par)
    Cap[i] = Retention_∂θ∂ψ(ψ[i], par)
  end
end

## Julia无法推测soil.param.param的类型
function cal_θKCap!(soil::Soil{T,P}, ψ::AbstractVector{T}) where {T<:Real,P<:AbstractSoilParam{T}}
  cal_θKCap!(soil, soil.param.param, ψ)
end

@inline function K_interface(K1::T, K2::T, d1::T, d2::T)::T where {T<:Real}
  # K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  return K1 * K2 * (d1 + d2) / (K1 * d2 + K2 * d1) # harmonic mean, # Eq. 5.16
end

function cal_K₊ₕ!(soil::Soil{T}) where {T<:Real}
  (; N, ibeg, Δz, K) = soil
  K₊ₕ::Vector{T} = soil.K₊ₕ
  @inbounds for i = ibeg:N-1
    K₊ₕ[i] = K_interface(K[i], K[i+1], Δz[i], Δz[i+1])
  end
end

cal_K!(soil::Soil) = cal_K!(soil, soil.θ)
# cal_K!(soil::Soil, θ::AbstractVector{T}) where {T<:Real} = cal_K!(soil, soil.param.param, θ)
function cal_K!(soil::Soil{T,P}, θ::AbstractVector{T}) where {T<:Real, P<:AbstractSoilParam{T}}
  (; N, ibeg, K) = soil
  param = soil.param.param
  @inbounds for i = ibeg:N
    K[i] = Retention_K(θ[i], param[i])
  end
end


cal_ψ!(soil::Soil) = cal_ψ!(soil, soil.θ)

function cal_ψ!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, ψ) = soil
  param = soil.param.param
  @inbounds for i = ibeg:N
    ψ[i] = Retention_ψ(θ[i], param[i])
  end
end


function Init_ψ0(soil::Soil{T}, θ::T) where {T<:Real}
  i = soil.ibeg # 默认采用第一层进行初始化
  param = soil.param.param
  return Retention_ψ(θ, param[i]) # ψ0
end



export specific_yield!
export cal_K!, cal_ψ!, cal_θKCap!, cal_K₊ₕ!
export Init_SoilWaterParam, Init_ψ0
