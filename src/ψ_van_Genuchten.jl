"""
    van_Genuchten(ψ, param)

van Genuchten (1980) relationships

# Arguments
+ `ψ`: Matric potential
+ `param`
  - `θ_res`       : Residual water content
  - `θ_sat`       : Volumetric water content at saturation
  - `α`           : Inverse of the air entry potential (cm-1)
  - `n`           : Pore-size distribution index
  - `m`           : Exponent
  - `K_sat`       : Hydraulic conductivity at saturation (cm/s)
  - `soil_texture`: Soil texture flag

# Examples
```julia
# Haverkamp et al. (1977): sand
param = (soil_texture = 1, 
  θ_res = 0.075, θ_sat = 0.287, 
  α = 0.027, n = 3.96, m = 1, K_sat = 34 / 3600)

# Haverkamp et al. (1977): Yolo light clay
param = (soil_texture=2, 
  θ_res = 0.124, θ_sat = 0.495,
  α = 0.026, n = 1.43, m = 1 - 1 / 1.43,
  K_sat = 0.0443 / 3600)
```
"""
function van_Genuchten(ψ::T; param::ParamVanGenuchten{T}, soil_type::Int=1) where {T<:Real}
  @unpack θ_res, θ_sat, α, n, m, Ksat = param

  # Effective saturation (Se) for specified matric potential (ψ)
  Se = ψ <= 0 ? (1 + (α * abs(ψ))^n)^-m : 1

  # Volumetric soil moisture (θ) for specified matric potential (ψ)
  θ = θ_res + (θ_sat - θ_res) * Se
  K = Ksat # 饱和时的水力传导率
  # Hydraulic conductivity (K) for specified matric potential (ψ)
  if Se <= 1
    K = Ksat * sqrt(Se) * (1 - (1 - Se^(1 / m))^m)^2
    # Special case for Haverkamp et al. (1977) sand (soil_texture = 1) and Yolo light clay (soil_texture = 2)
    if soil_type == 1
      K = Ksat * 1.175e6 / (1.175e6 + abs(ψ)^4.74)
    elseif soil_type == 2
      K = Ksat * 124.6 / (124.6 + abs(ψ)^1.77)
    end
  end

  # Specific moisture capacity (∂θ∂ψ) for specified matric potential (ψ)
  if ψ <= 0.0
    num = α * m * n * (θ_sat - θ_res) * (α * abs(ψ))^(n - 1)
    den = (1 + (α * abs(ψ))^n)^(m + 1)
    ∂θ∂ψ = num / den
  else
    ∂θ∂ψ = 0.0
  end
  θ, K, ∂θ∂ψ
end


# Function to calculate hydraulic conductivity from water content
function van_genuchten_K(θ::T; param::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, Ksat) = param
  Se = (θ - θ_res) / (θ_sat - θ_res)
  # Se = clamp(Se, 0, 1)
  # effective_saturation = Se^0.5
  # term = (1 - (1 - Se^(1 / m))^m)^2
  # return Ksat * effective_saturation * term

  if Se <= 1
    # Special case for:
    # - `soil_texture = 1`: Haverkamp et al. (1977) sand
    # - `soil_texture = 2`: Yolo light clay
    ψ = van_genuchten_ψ(θ; param)
    return Ksat * 1.175e6 / (1.175e6 + pow(abs(ψ), 4.74)) # Haverkamp et al. (1977) sand
  # Ksat * 124.6 / (124.6 + abs(ψ)^1.77)   # Yolo light clay
  else
    return Ksat
  end
end

# Function to calculate pressure head psi from water content
function van_genuchten_ψ(θ::T; param::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, α, n, m) = param
  if θ <= θ_res
    return T(-Inf)  # Return a very high positive number indicating very dry conditions
  elseif θ >= θ_sat
    return T(0.0)  # Saturated condition, psi is zero
  else
    return -1 / α * pow(pow((θ_sat - θ_res) / (θ - θ_res), (1 / m)) - 1, 1 / n)
  end
end
