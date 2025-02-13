"""
  计算每一层的给水度

Sy = θ_sat - θ(ψ_sat + zwt)

前提条件是：每层土壤均接近饱和，才能得到这个公式，`ψ + z = ψ_sat + zwt`，即Q = 0。
"""
function specific_yield!(soil::Soil{T}, zwt::T; sy_max::T=0.02) where {T<:Real}
  (; N, Sy) = soil
  (; θ_sat, ψ_sat, method) = soil.param
  (; param) = soil.param

  iszero_ψ = method == "van_Genuchten"
  for i = 1:N
    _ψ_sat = iszero_ψ ? T(0) : ψ_sat[i]
    Sy[i] = θ_sat[i] - Retention_θ(_ψ_sat + zwt, param[i])
  end
  clamp!(Sy, 0, sy_max)
end


function cal_θKCap!(soil::Soil{T}, ψ::AbstractVector{T}) where {T<:Real}
  (; ibeg, N, θ, K, Cap) = soil
  (; param) = soil.param
  @inbounds for i in ibeg:N
    θ[i], K[i], Cap[i] = Retention(ψ[i], param[i])
  end
end


function cal_K₊ₕ!(soil::Soil)
  (; N, ibeg, z, z₊ₕ, K, K₊ₕ) = soil
  @inbounds for i = ibeg:N-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    K₊ₕ[i] = K[i] * K[i+1] * (d1 + d2) / (K[i] * d2 + K[i+1] * d1) # Eq. 5.16, 
    # K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
end

cal_K!(soil::Soil) = cal_K!(soil, soil.θ)
function cal_K!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, K) = soil
  (; param) = soil.param
  @inbounds for i = ibeg:N
    K[i] = Retention_K(θ[i], param[i])
  end
end


cal_ψ!(soil::Soil) = cal_ψ!(soil, soil.θ)
function cal_ψ!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, ψ) = soil
  (; param) = soil.param
  @inbounds for i = ibeg:N
    ψ[i] = Retention_ψ(θ[i], param[i])
  end
end


function Init_ψ0(soil::Soil{T}, θ::T) where {T<:Real}
  i = soil.ibeg # 默认采用第一层进行初始化
  (; param) = soil.param
  return Retention_ψ(θ, param[i]) # ψ0
end

function Init_SoilWaterParam(N, θ_sat::T, θ_res::T, Ksat::T, α::T, n::T, m::T; use_m::Bool=false, same_layer=true, kw...) where {T<:Real}
  SoilParam(; N,
    θ_sat=fill(θ_sat, N),
    θ_res=fill(θ_res, N),
    Ksat=fill(Ksat, N),
    α=fill(α, N),
    n=fill(n, N),
    m=fill(m, N),
    use_m, same_layer, kw...)
end


export specific_yield!
export cal_K!, cal_ψ!, cal_θKCap!, cal_K₊ₕ!
export Init_SoilWaterParam, Init_ψ0
