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
  K::Vector{FT} = zeros(FT, n)       # hydraulic conductivity，[cm/s]
  K₊ₕ::Vector{FT} = zeros(FT, n - 1)  # hydraulic conductivity at interface, [cm/s]
  Cap::Vector{FT} = zeros(FT, n)     # specific moisture capacity, dθ/dΨ, [cm-1], 临时变量
  ψ::Vector{FT} = zeros(FT, n)       # [cm]，约干越负
  ψ_next::Vector{FT} = zeros(FT, n)  # ψ[n+1/2], [cm], 临时变量
  θ0::FT = FT(0.0)                   # [m3 m-3]
  ψ0::FT = FT(0.0)                   # [cm]
  Q0::FT = FT(0.0)                   # [cm/s] 下渗速率，向下为负
  sink::Vector{FT} = fill(0.0, n)    # 蒸发项, [cm per unit time]
  θ_prev::Vector{FT} = zeros(FT, n)  # backup of θ
  ψ_prev::Vector{FT} = zeros(FT, n)  # backup of ψ

  # 温度
  Tsoil::Vector{FT} = fill(NaN, n)   # [°C]
  κ::Vector{FT} = zeros(FT, n)       # thermal conductivity [W m-1 K-1]
  κ₊ₕ::Vector{FT} = zeros(FT, n - 1)  # thermal conductivity at interface [W m-1 K-1]
  cv::Vector{FT} = zeros(FT, n)      # volumetric heat capacity [J m-3 K-1]
  F::Vector{FT} = zeros(FT, n)       # heat flux, [W m-2]
  TS0::FT = FT(NaN)                  # surface temperature, [°C]
  F0::FT = FT(NaN)                   # heat flux at the surface, [W m-2]，向下为负
  G::FT = FT(NaN)                    # [W m-2]，土壤热通量

  # ODE求解临时变量
  u::Vector{FT} = fill(NaN, n)  # [°C], 为了从ibeg求解地温，定义的临时变量
  du::Vector{FT} = fill(NaN, n) # [°C]

  # 三角阵求解临时变量
  a::Vector{FT} = zeros(FT, n)
  b::Vector{FT} = zeros(FT, n)
  c::Vector{FT} = zeros(FT, n)
  d::Vector{FT} = zeros(FT, n)
  e::Vector{FT} = zeros(FT, n)
  f::Vector{FT} = zeros(FT, n)

  timestep::Int = 0                  # 迭代次数
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
