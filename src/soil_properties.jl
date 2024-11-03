
# update θ, K, Cap
function cal_θKCap!(soil::Soil{T}, ψ::AbstractVector{T}) where {T<:Real}
  (; ibeg, N, θ, K, Cap) = soil
  (; θ_sat, θ_res, Ksat, α, n, m, use_m, b, method) = soil.param

  if method == "van_Genuchten"
    @inbounds for i in ibeg:N
      # _m = m[i]
      _m = use_m ? m[i] : 1 - 1 / n[i] # 不参与参数优化
      θ[i], K[i], Cap[i] = van_Genuchten(ψ[i], θ_sat[i], θ_res[i], Ksat[i], α[i], n[i], _m)
    end
  elseif method == "Campbell"
    @inbounds for i in ibeg:N
      θ[i], K[i], Cap[i] = Cambell(ψ[i], ψ_sat[i], θ_sat[i], Ksat[i], b[i])
    end
  end
end


# update K₊ₕ
function cal_K₊ₕ!(soil::Soil)
  (; N, ibeg, z, z₊ₕ, K, K₊ₕ) = soil
  @inbounds for i = ibeg:N-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    K₊ₕ[i] = K[i] * K[i+1] * (d1 + d2) / (K[i] * d2 + K[i+1] * d1) # Eq. 5.16, 
    # K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
end


function cal_K!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, K) = soil
  (; θ_sat, θ_res, Ksat, m, b, method) = soil.param

  if method == "van_Genuchten"
    @inbounds for i = ibeg:N
      K[i] = van_Genuchten_K(θ[i], θ_sat[i], θ_res[i], Ksat[i], m[i])
    end
  elseif method == "Campbell"
    @inbounds for i = ibeg:N
      K[i] = Cambell_K(θ[i], θ_sat[i], Ksat[i], b[i])
    end
  end
end


# cal_ψ!(soil, θ)
function cal_ψ!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, ψ) = soil
  (; θ_sat, θ_res, α, n, m, ψ_sat, b, method) = soil.param
  if method == "van_Genuchten"
    @inbounds for i = ibeg:N
      ψ[i] = van_Genuchten_ψ(θ[i], θ_sat[i], θ_res[i], α[i], n[i], m[i])
    end
  elseif method == "Campbell"
    @inbounds for i = ibeg:N
      ψ[i] = Campbell_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
    end
  end
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

function Init_ψ0(soil::Soil{T}, θ::T) where {T<:Real}
  i = soil.ibeg # 默认采用第一层进行初始化
  (; θ_sat, θ_res, α, n, m, ψ_sat, b, method) = soil.param

  if method == "van_Genuchten"
    ψ0 = van_Genuchten_ψ(θ, θ_sat[i], θ_res[i], α[i], n[i], m[i])
  elseif method == "Campbell"
    ψ0 = Campbell_ψ(θ, θ_sat[i], ψ_sat[i], b[i])
  end
  return ψ0
end


export Init_SoilWaterParam, Init_ψ0
