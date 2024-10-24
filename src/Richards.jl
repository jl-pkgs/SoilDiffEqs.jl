using DifferentialEquations
using HydroTools
using Plots
using Parameters

include("Soil.jl")
include("soil_moisture_Q0.jl")
include("soil_moisture.jl")

"""
# Arguments

- `method`: 
  + `TS0`: TS0 boundary condition, 第一类边界条件
  + `F0`: F0 boundary condition, 第二类边界条件
"""
function SoilTemperature(dT, T, p::Soil, t; method="F0")
  (; TS0, κ, cv, z) = p
  if method == "TS0"
    F0 = -κ[1] * (TS0 - T[1])
  elseif method == "F0"
    F0 = p.F0
  end

  @inbounds for i in 1:n-1
    κ₊ₕ = (κ[i] + κ[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    F[i] = -κ₊ₕ * (T[i] - T[i+1]) / Δz₊ₕ
  end
  F[n] = 0

  dT[1] = -(F0 - F[1]) / (Δz[1] * cv[1])
  @inbounds for i in 2:n
    dT[i] = -(F[i-1] - F[i]) / (Δz[i] * cv[i])
  end
end

"""
# Arguments
- `method`: 
  + `ψ0`: ψ0 boundary condition, 第一类边界条件
  + `Q0`: Q0 boundary condition, 第二类边界条件
"""
function RichardsEquation(dθ, u, p::Soil, t; method="ψ0")
  (; n, Q, z, Δz, K, ψ, ψ0, sink) = p
  param = p.param_water
  @. K = van_genuchten_K(u; param)
  @. ψ = van_genuchten_ψ(u; param)

  if method == "ψ0"
    Q0 = -K[1] * ((ψ0 - ψ[1]) / (0 - z[1]) + 1)
  elseif method == "Q0"
    Q0 = p.Q0
  end

  @inbounds for i in 1:n-1
    K₊ₕ = (K[i] + K[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    Q[i] = -K₊ₕ * ((ψ[i] - ψ[i+1]) / Δz₊ₕ + 1)
  end
  Q[n] = -K[n]

  dθ[1] = -(Q0 - Q[1]) / Δz[1] - sink[1] / Δz[1]
  @inbounds for i in 2:n
    dθ[i] = -(Q[i-1] - Q[i]) / Δz[i] - sink[i] / Δz[i]
  end
end

# Q[i-1]
function cal_Q(i::Int, K::Vector{FT}, ψ::Vector{FT}, z::Vector{FT}) where {FT}
  K₊ₕ = (K[i] + K[i-1]) / 2
  Δz₊ₕ = z[i-1] - z[i]
  Q = -K₊ₕ * ((ψ[i-1] - ψ[i]) / Δz₊ₕ + 1)
  Q
end

# function RichardsEquation_low(du, u, p::Soil, t)
#   (; n, z, z₊ₕ, K, ψ, ψ0, param) = p

#   @. K = van_genuchten_K(u; param)
#   @. ψ = van_genuchten_ψ(u; param)

#   # @inbounds 
#   @inbounds for i in 2:n-1
#     Q_up = cal_Q(i, K, ψ, z) # 多计算了1次
#     Q_down = cal_Q(i + 1, K, ψ, z)
#     Δz = z₊ₕ[i-1] - z₊ₕ[i] # 
#     du[i] = -(Q_up - Q_down) / Δz
#   end

#   ## boundary
#   i = 1
#   K₊ₕ_up = K[1]
#   Q_up = -K₊ₕ_up * ((ψ0 - ψ[i]) / (0 - z[i]) + 1)
#   Q_down = cal_Q(i + 1, K, ψ, z)
#   Δz = 0 - z₊ₕ[i]
#   du[i] = -(Q_up - Q_down) / Δz

#   i = n
#   Q_up = cal_Q(i, K, ψ, z)
#   Q_down = -K[i]
#   Δz = (z[i-1] - z[i])
#   du[i] = -(Q_up - Q_down) / Δz
# end
