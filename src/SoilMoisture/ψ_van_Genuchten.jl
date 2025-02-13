"""
    Retention(ψ::T, par::ParamVanGenuchten{T})

van Genuchten (1980) relationships

# Arguments
+ `ψ`: Matric potential
+ `par`
  - `θ_res`       : Residual water content
  - `θ_sat`       : Volumetric water content at saturation
  - `Ksat`        : Hydraulic conductivity at saturation (cm/s)
  - `α`           : Inverse of the air entry potential (cm-1)
  - `n`           : Pore-size distribution index
  - `m`           : Exponent
  - `soil_texture`: Soil texture flag

# Examples
```julia
# Haverkamp et al. (1977): sand
par = (soil_texture = 1, 
  θ_res = 0.075, θ_sat = 0.287, 
  α = 0.027, n = 3.96, m = 1, Ksat = 34 / 3600)

# Haverkamp et al. (1977): Yolo light clay
par = (soil_texture=2, 
  θ_res = 0.124, θ_sat = 0.495,
  α = 0.026, n = 1.43, m = 1 - 1 / 1.43,
  Ksat = 0.0443 / 3600)
```
"""
function Retention_VanGenuchten(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_res, θ_sat, Ksat, α, n, m) = par
  # θ_sat::T, θ_res::T, Ksat::T, α::T, n::T, m::T
  # Effective saturation (Se) for specified matric potential (ψ)
  Se = ψ <= 0 ? (1 + (α * abs(ψ))^n)^-m : 1

  # Volumetric soil moisture (θ) for specified matric potential (ψ)
  θ = θ_res + (θ_sat - θ_res) * Se
  # Hydraulic conductivity (K) for specified matric potential (ψ)
  diff = (1.0 - Se^(1 / m))
  K = Se < 1 ? Ksat * sqrt(Se) * (1 - diff^m)^2 : Ksat

  # Specific moisture capacity (∂θ∂ψ) for specified matric potential (ψ)
  ∂θ∂ψ = Retention_∂θ∂ψ(ψ, par)
  θ, K, ∂θ∂ψ
end

@inline @fastmath function Retention_∂θ∂ψ(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_res, θ_sat, α, n, m) = par
  if ψ <= 0.0
    num = α * m * n * (θ_sat - θ_res) * (α * abs(ψ))^(n - 1)
    den = (1 + (α * abs(ψ))^n)^(m + 1)
    ∂θ∂ψ = num / den
  else
    ∂θ∂ψ = 0.0
  end
end

function Retention_θ(ψ::T, par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, α, n, m) = par
  # Effective saturation (Se) for specified matric potential (ψ)
  Se = ψ <= 0 ? (1 + (α * abs(ψ))^n)^-m : 1
  # Volumetric soil moisture (θ) for specified matric potential (ψ)
  return θ_res + (θ_sat - θ_res) * Se # θ
end

function Retention_K(θ::T, par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, Ksat, m) = par
  Se = (θ - θ_res) / (θ_sat - θ_res)
  Se = clamp(Se, 0.0, 1.0)
  diff = (1.0 - Se^(1 / m))
  K = Se < 1 ? Ksat * sqrt(Se) * (1 - diff^m)^2 : Ksat
  return K
end

function Retention_ψ(θ::T, par::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, α, n, m) = par
  if θ <= θ_res
    return T(-Inf)  # Return a very high positive number indicating very dry conditions
  elseif θ >= θ_sat
    return T(0.0)   # Saturated condition, psi is zero
  else
    return -1 / α * pow(pow((θ_sat - θ_res) / (θ - θ_res), (1 / m)) - 1, 1 / n)
  end
end

Retention_ψ(θ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_ψ(θ, par)
Retention_K(θ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_K(θ, par)
Retention_θ(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_θ(ψ, par)

# Special case for:
# - `soil_texture = 1`: Haverkamp et al. (1977) sand
# - `soil_texture = 2`: Yolo light clay

# if soil_type == 1
#   K = Ksat * 1.175e6 / (1.175e6 + abs(ψ)^4.74)
# elseif soil_type == 2
#   K = Ksat * 124.6 / (124.6 + abs(ψ)^1.77)
# end
