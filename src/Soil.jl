using Parameters

# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = x^y
# @fastmath pow(x::Real, y::Real) = exp(y * log(x))

# Ksat: [cm/s]
abstract type AbstractSoilParam{FT} end

@with_kw mutable struct ParamVanGenuchten{T} <: AbstractSoilParam{T}
  θs::T = 0.287       # [m3 m-3]
  θr::T = 0.075       # [m3 m-3]
  Ksat::T = 34 / 3600 # [cm s-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0
end


# 两种边界条件的设定方法：
@with_kw mutable struct Soil{FT}
  n::Int = 10                        # layers of soil
  dt::Float64 = 3600                 # 时间步长, seconds
  z::Vector{FT} = zeros(FT, n)       # cm, 向下为负
  z₊ₕ::Vector{FT} = zeros(FT, n)
  Δz::Vector{FT} = zeros(FT, n)
  Δz₊ₕ::Vector{FT} = zeros(FT, n)

  # 水分
  θ::Vector{FT} = fill(0.1, n)       # θ [m3 m-3]
  Q::Vector{FT} = zeros(FT, n)       # [cm/s]
  K::Vector{FT} = zeros(FT, n)       # 水力传导系数，[cm/s]
  ψ::Vector{FT} = zeros(FT, n)       # [cm]，约干越负
  ψ0::FT = FT(0.0)                   # [cm]
  Q0::FT = FT(0.0)                   # [cm/s] 下渗速率，向下为负
  sink::Vector{FT} = fill(0.0, n)    # 蒸发项, [cm per unit time]

  # 温度
  Tsoil::Vector{FT} = ones(FT, n) .* NaN # [°C]
  κ::Vector{FT} = zeros(FT, n)       # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = zeros(FT, n)      # volumetric heat capacity [J m-3 K-1]
  F::Vector{FT} = zeros(FT, n)       # heat flux, [W m-2]
  TS0::FT = FT(NaN)                  # surface temperature, [°C]
  F0::FT = FT(NaN)                   # heat flux at the surface, [W m-2]

  timestep::Int = 0             # 时间步长
  param_water::ParamVanGenuchten{FT} = ParamVanGenuchten{FT}()
end

# Function to calculate hydraulic conductivity from water content
function van_genuchten_K(θ::T; param::ParamVanGenuchten{T}) where {T<:Real}
  (; θs, θr, Ksat) = param
  Se = (θ - θr) / (θs - θr)
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
  (; θs, θr, α, n, m) = param
  if θ <= θr
    return -Inf  # Return a very high positive number indicating very dry conditions
  elseif θ >= θs
    return 0  # Saturated condition, psi is zero
  else
    return -1 / α * pow(pow((θs - θr) / (θ - θr), (1 / m)) - 1, 1 / n)
  end
end
