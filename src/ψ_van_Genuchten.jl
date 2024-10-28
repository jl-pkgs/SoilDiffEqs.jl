export get_soilpar, ParamVanGenuchten

# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = x^y
# @fastmath pow(x::Real, y::Real) = exp(y * log(x))

# Ksat: [cm/s]
abstract type AbstractSoilParam{FT} end

@with_kw mutable struct ParamVanGenuchten{T} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  θ_res::T = 0.075       # [m3 m-3]
  Ksat::T = 34 / 3600 # [cm s-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0 - 1.0 / n
end

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
  Ksat = 0.0443 / 3600)
```
"""
function van_Genuchten(ψ::T; param::ParamVanGenuchten{T}) where {T<:Real}
  @unpack θ_res, θ_sat, α, n, m, Ksat = param

  # Effective saturation (Se) for specified matric potential (ψ)
  Se = ψ <= 0 ? (1 + (α * abs(ψ))^n)^-m : 1

  # Volumetric soil moisture (θ) for specified matric potential (ψ)
  θ = θ_res + (θ_sat - θ_res) * Se

  # Hydraulic conductivity (K) for specified matric potential (ψ)
  K = Se < 1 ? Ksat * sqrt(Se) * (1 - (1 - Se^(1 / m))^m)^2 : Ksat

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
  (; θ_sat, θ_res, Ksat, m) = param
  Se = (θ - θ_res) / (θ_sat - θ_res)
  K = Se < 1 ? Ksat * sqrt(Se) * (1 - (1 - Se^(1 / m))^m)^2 : Ksat
  return K
end

# Function to calculate pressure head psi from water content
function van_genuchten_ψ(θ::T; param::ParamVanGenuchten{T}) where {T<:Real}
  (; θ_sat, θ_res, α, n, m) = param
  if θ <= θ_res
    return T(-Inf)  # Return a very high positive number indicating very dry conditions
  elseif θ >= θ_sat
    return T(0.0)   # Saturated condition, psi is zero
  else
    return -1 / α * pow(pow((θ_sat - θ_res) / (θ - θ_res), (1 / m)) - 1, 1 / n)
  end
end

# Special case for:
# - `soil_texture = 1`: Haverkamp et al. (1977) sand
# - `soil_texture = 2`: Yolo light clay

# if soil_type == 1
#   K = Ksat * 1.175e6 / (1.175e6 + abs(ψ)^4.74)
# elseif soil_type == 2
#   K = Ksat * 124.6 / (124.6 + abs(ψ)^1.77)
# end


function Base.Vector(x::ParamVanGenuchten)
  (; θ_sat, θ_res, Ksat, α, n, m) = x
  [θ_sat, θ_res, Ksat, α, n, m]
end

# Bonan 2019, Table 8.3
function get_soilpar(soil_type::Int=1)
  soilparam = [
    # θ_sat, θ_res, α (cm⁻¹), n, Ksat (cm h⁻¹)
    0.38 0.068 0.008 1.09 0.2;   #  1,  Clay
    0.36 0.07 0.005 1.09 0.02;  #  2,  Silty  clay
    0.38 0.1 0.027 1.23 0.12;  #  3,  Sandy  clay
    0.41 0.095 0.019 1.31 0.26;  #  4,  Clay   loam
    0.43 0.089 0.01 1.23 0.07;  #  5,  Silty  clay loam
    0.39 0.1 0.059 1.48 1.31;  #  6,  Sandy  clay loam
    0.43 0.078 0.036 1.56 1.04;  #  7,  Loam
    0.45 0.067 0.02 1.41 0.45;  #  8,  Silty  loam
    0.41 0.065 0.075 1.89 4.42;  #  9,  Sandy  loam
    0.41 0.065 0.075 1.89 4.42;  #  10, Silty, no   data in Bonan2019
    0.41 0.057 0.124 2.28 14.59; #  11, Loamy  sand
    0.43 0.045 0.145 2.68 29.7   #  12, Sand
  ]
  θ_sat, θ_res, α, n, Ksat = soilparam[soil_type, :]
  Ksat = Ksat / 3600 # [cm h-1] to [cm s-1]
  ParamVanGenuchten(; θ_sat, θ_res, α, n, Ksat)
end

function get_soilpar(theta::AbstractVector)
  θ_sat, θ_res, Ksat, α, n = theta[1:5]
  ParamVanGenuchten(; θ_sat, θ_res, α, n, Ksat)
end
