"""
    Cambell(ψ, ψ_sat, θ_sat, Ksat, b)

Campbell (1974) relationships

# Arguments
+ `ψ`: Matric potential, cm
+ `param`
  - `ψ_sat`: Matric potential at saturation, [cm]
  - `θ_sat`: Volumetric water content at saturation
  - `b`    : Exponent
  - `Ksat`: Hydraulic conductivity at saturation [cm h-1]

# TODO: 核对变量的单位

# Examples
```julia
param = (θ_sat = 0.25, ψ_sat = -25.0, b = 0.2, Ksat = 3.4e-03)
θ, K, ∂θ∂ψ = Cambell(ψ, ψ_sat, θ_sat, Ksat, b)
```
"""
@fastmath function Cambell(ψ::T, ψ_sat::T, θ_sat::T, Ksat::T, b::T) where {T<:Real}
  # @unpack ψ_sat, θ_sat, Ksat, b = param

  # Volumetric soil moisture (θ) for specified matric potential 
  θ = ψ < ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : ψ_sat

  # Hydraulic conductivity (K) for specified matric potential
  K = ψ < ψ_sat ? Ksat * (θ / θ_sat)^(2b + 3) : Ksat

  # Cap = dθ/dψ
  ∂θ∂ψ = ψ < ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : ψ_sat

  θ, K, ∂θ∂ψ
end

"""
    Cambell_θ(ψ, ψ_sat, θ_sat, b)
"""
@inline function Cambell_θ(ψ::T, ψ_sat::T, θ_sat::T, b::T) where {T<:Real}
  return ψ < ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : ψ_sat # θ
end

"""
    Cambell_K(θ, θ_sat, Ksat, b)
"""
@fastmath Cambell_K(θ::T, θ_sat::T, Ksat::T, b::T) where {T<:Real} =
  Ksat * (θ / θ_sat)^(2 * b + 3)

# Campbell 1974, Bonan 2019 Table 8.2
"""
    Cambell_ψ(θ, θ_sat, ψ_sat, b)
"""
@fastmath function Cambell_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end
