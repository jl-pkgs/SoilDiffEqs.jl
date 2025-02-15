export Retention, Retention_K, Retention_θ, Retention_ψ, Retention_∂θ∂ψ, Retention_∂K∂Se

# 多重派发(runtime-dispatch)可能会导致速度变慢
Retention(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell(ψ, par)
Retention_K(θ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_K(θ, par)
Retention_θ(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_θ(ψ, par)
Retention_ψ(θ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_ψ(θ, par)
Retention_ψ_Se(Se::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_ψ_Se(Se, par)
Retention_∂θ∂ψ(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_∂θ∂ψ(ψ, par)
Retention_∂ψ∂θ(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_∂ψ∂θ(ψ, par)
Retention_∂K∂Se(Se::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_∂K∂Se(Se, par)


Retention(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten(ψ, par)
Retention_θ(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_θ(ψ, par)
Retention_K(θ::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_K(θ, par)
Retention_ψ(θ::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_ψ(θ, par)
Retention_ψ_Se(Se::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_ψ_Se(Se, par)
Retention_∂θ∂ψ(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_∂θ∂ψ(ψ, par)
Retention_∂ψ∂θ(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_∂ψ∂θ(ψ, par)
Retention_∂K∂Se(Se::T, par::ParamVanGenuchten{T}) where {T<:Real} = van_Genuchten_∂K∂Se(Se, par)


Retention(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention(ψ, par)
Retention_K(θ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_K(θ, par)
Retention_θ(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_θ(ψ, par)
Retention_ψ(θ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_ψ(θ, par)
Retention_ψ_Se(Se::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_ψ_Se(Se, par)
Retention_∂θ∂ψ(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_∂θ∂ψ(ψ, par)
Retention_∂ψ∂θ(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_∂ψ∂θ(ψ, par)
Retention_∂K∂Se(Se::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_∂K∂Se(Se, par)


"""
  计算每一层的给水度

Sy = θ_sat - θ(ψ_sat + zwt)

前提条件是：每层土壤均接近饱和，才能得到这个公式，`ψ + z = ψ_sat + zwt`，即Q = 0。
"""
function specific_yield!(soil::Soil{T}, zwt::T; sy_max::T=0.02) where {T<:Real}
  (; N, Sy) = soil
  (; θ_sat, param, method_retention) = soil.param

  iszero_ψ = method_retention == "van_Genuchten"
  for i = 1:N
    _ψ_sat = iszero_ψ ? 0.0 : θ_sat[i] # !Note
    Sy[i] = θ_sat[i] - Retention_θ(_ψ_sat + zwt, param[i])
  end
  clamp!(Sy, 0.0, sy_max)
end


function cal_θ!(soil::Soil, ψ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, θ) = soil
  (; param) = soil.param
  @inbounds for i = ibeg:N
    θ[i] = Retention_θ(ψ[i], param[i])
  end
end

function cal_∂θ∂ψ!(soil::Soil{T,P}, ψ::AbstractVector{T}) where {T<:Real,P<:AbstractSoilParam{T}}
  (; ibeg, N, ∂θ∂ψ) = soil
  (; param) = soil.param
  @inbounds for i in ibeg:N
    ∂θ∂ψ[i] = Retention_∂θ∂ψ(ψ[i], param[i])
  end
end

# 算术平均、调和平均
mean_arithmetic(K1::T, K2::T, d1::T, d2::T) where {T<:Real} = (K1 * d1 + K2 * d2) / (d1 + d2)
mean_harmonic(K1::T, K2::T, d1::T, d2::T) where {T<:Real} = K1 * K2 * (d1 + d2) / (K1 * d2 + K2 * d1)

# 一并更新K和K₊ₕ
cal_K!(soil::Soil) = cal_K!(soil, soil.θ)
function cal_K!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, K, K₊ₕ, Δz) = soil
  param = soil.param.param
  i0 = max(soil.ibeg - 1, 1)
  
  @inbounds for i = i0:N
    K[i] = Retention_K(θ[i], param[i])
  end
  @inbounds for i = i0:N-1
    K₊ₕ[i] = mean_arithmetic(K[i], K[i+1], Δz[i], Δz[i+1])
    # K₊ₕ[i] = mean_harmonic(K[i], K[i+1], Δz[i], Δz[i+1])
  end
  K₊ₕ[N] = K[N]
end


cal_ψ!(soil::Soil) = cal_ψ!(soil, soil.θ)
function cal_ψ!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ψ) = soil
  param = soil.param.param
  i0 = max(soil.ibeg - 1, 1)
  
  @inbounds for i = i0:N
    ψ[i] = Retention_ψ(θ[i], param[i])
  end
end


function Init_ψ0(soil::Soil{T}, θ::T) where {T<:Real}
  # ibeg = 1, ψ0; ibeg = 2, ψ1
  i = max(soil.ibeg - 1, 1)
  return Retention_ψ(θ, soil.param.param[i]) # ψ0
end


export specific_yield!
export cal_K!, cal_ψ!, cal_θKCap!, cal_K₊ₕ!
export Init_SoilWaterParam, Init_ψ0
