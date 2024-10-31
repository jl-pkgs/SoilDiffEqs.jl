# update θ, K, Cap
function update_θ!(soil::Soil{T}, ψ::AbstractVector{T}) where {T<:Real}
  (; ibeg, N, θ, K, Cap) = soil
  (; θ_sat, θ_res, Ksat, α, n, b, method) = soil.param

  if method == "van_Genuchten"
    for i in ibeg:N
      m = 1 - 1 / n[i] # m，不参与参数优化
      θ[i], K[i], Cap[i] = van_Genuchten(ψ[i], θ_res[i], θ_sat[i], Ksat[i], α[i], n[i], m)
      # θ[i], K[i], Cap[i] = van_Genuchten(ψ[i]; param)
    end
  elseif method == "Campbell"
    for i in ibeg:N
      θ[i], K[i], Cap[i] = Cambell(ψ[i], ψ_sat[i], θ_sat[i], Ksat[i], b[i])
    end
  end
end


# update K₊ₕ
function update_K₊ₕ!(soil::Soil)
  (; N, ibeg, z, z₊ₕ, K, K₊ₕ) = soil
  for i = ibeg:N-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    K₊ₕ[i] = K[i] * K[i+1] * (d1 + d2) / (K[i] * d2 + K[i+1] * d1) # Eq. 5.16, 
    # K₊ₕ[i] = (K[i] + K[i+1]) / 2 # can be improved, weighted by z
  end
end


function update_K!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, K) = soil
  (; θ_sat, θ_res, Ksat, m, b, method) = soil.param

  if method == "van_Genuchten"
    for i = ibeg:N
      K[i] = van_Genuchten_K(θ[i], θ_sat[i], θ_res[i], Ksat[i], m[i])
    end
  elseif method == "Campbell"
    for i = ibeg:N
      K[i] = Cambell_K(θ[i], θ_sat[i], Ksat[i], b[i])
    end
  end
end


function update_ψ!(soil::Soil, θ::AbstractVector{T}) where {T<:Real}
  (; N, ibeg, ψ) = soil
  (; θ_sat, θ_res, α, n, m, b, method) = soil.param
  if method == "van_Genuchten"
    for i = ibeg:N
      ψ[i] = van_Genuchten_ψ(θ[i], θ_sat[i], θ_res[i], α[i], n[i], m[i])
    end
  elseif method == "Campbell"
    for i = ibeg:N
      ψ[i] = Campbell_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
    end
  end
end
