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
function SoilTemperature(dT, T, p::Soil, t; method="TS0")
  p.timestep += 1
  (; n, Δz, z, z₊ₕ, F, κ, cv, TS0) = p

  if method == "TS0"
    F0 = -κ[1] * (TS0 - T[1]) / (0 - z[1])
  elseif method == "F0"
    F0 = p.F0
  end
  
  @inbounds for i in 1:n-1
    κ₊ₕ = κ[i] * κ[i+1] * (z[i] - z[i+1]) /
             (κ[i] * (z₊ₕ[i] - z[i+1]) + κ[i+1] * (z[i] - z₊ₕ[i])) # Eq. 5.16, 
    # κ₊ₕ = (κ[i] + κ[i+1]) / 2
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
function RichardsEquation(dθ::AbstractVector{T}, u::AbstractVector{T}, p::Soil{T}, t; method="ψ0") where {T<:Real}
  p.timestep += 1
  (; n, Δz, z, Q, K, ψ, ψ0, sink) = p
  param = p.param_water
  @inbounds for i = 1:n
    K[i] = van_genuchten_K(u[i]; param)
    ψ[i] = van_genuchten_ψ(u[i]; param)
  end
  # @. K = van_genuchten_K(u; param)
  # @. ψ = van_genuchten_ψ(u; param)
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
