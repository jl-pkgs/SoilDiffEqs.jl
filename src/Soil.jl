using DifferentialEquations
using HydroTools
using Plots
using Parameters
include("Soil_depth.jl")


# 两种边界条件的设定方法：
# 1. Q0
# 2. ψ0
@with_kw mutable struct Soil{FT}
  n::Int = 10
  z::Vector{FT} = zeros(FT, n)       # cm, 向下为负
  z₊ₕ::Vector{FT} = zeros(FT, n)
  Δz::Vector{FT} = zeros(FT, n)
  Δz₊ₕ::Vector{FT} = zeros(FT, n)

  # 水分
  θ::Vector{FT} = ones(FT, n) .* 0.1 # θ [m3 m-3]
  Q::Vector{FT} = zeros(FT, n)       # [cm/s]
  K::Vector{FT} = zeros(FT, n)       # 水力传导系数，[cm/s]
  ψ::Vector{FT} = zeros(FT, n)       # [cm]，约干越负
  ψ0::FT = FT(0.0)                   # [cm]
  Q0::FT = FT(0.0)                   # [cm/s] 下渗速率，向下为负
  sink::Vector{FT} = ones(FT, n) .* 0.0    # 蒸发项, [cm per unit time]

  # 温度
  κ::Vector{FT} = zeros(FT, n)
  cv::Vector{FT} = zeros(FT, n)
  TS0::FT = FT(NaN)                  # 边界层条件，地表温度
  F0::FT = FT(NaN)                   # 边界层条件，地表热通量

  param_water::NamedTuple = (; θs=0.287, θr=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1)
end
# Ksat: [cm/s]

# Function to calculate hydraulic conductivity from water content
function van_genuchten_K(θ; param)
  (; θs, θr, Ksat, α, n, m) = param
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
    return Ksat * 1.175e6 / (1.175e6 + abs(ψ)^4.74) # Haverkamp et al. (1977) sand
  # Ksat * 124.6 / (124.6 + abs(ψ)^1.77)   # Yolo light clay
  else
    return Ksat
  end
end

# Function to calculate pressure head psi from water content
function van_genuchten_ψ(θ; param)
  (; θs, θr, α, n, m) = param
  if θ <= θr
    return -Inf  # Return a very high positive number indicating very dry conditions
  elseif θ >= θs
    return 0  # Saturated condition, psi is zero
  else
    return -(1 / α) * (((θs - θr) / (θ - θr))^(1 / m) - 1)^(1 / n)
  end
end
