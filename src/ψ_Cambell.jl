"""
    Cambell(ψ; param)

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
θ, K, ∂θ∂ψ = Cambell(-100; param)
```
"""
@fastmath function Cambell(ψ::Real; param)
  @unpack ψ_sat, θ_sat, b, Ksat = param

  # Volumetric soil moisture (θ) for specified matric potential 
  θ = ψ < ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : ψ_sat

  # Hydraulic conductivity (K) for specified matric potential
  K = ψ < ψ_sat ? Ksat * (θ / θ_sat)^(2b + 3) : Ksat

  # Cap = dθ/dψ
  ∂θ∂ψ = ψ < ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : ψ_sat

  θ, K, ∂θ∂ψ
end

# Campbell 1974, Bonan 2019 Table 8.2
@fastmath function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end

@fastmath cal_K(θ::T, θ_sat::T, Ksat::T, b::T) where {T<:Real} =
  Ksat * (θ / θ_sat)^(2 * b + 3)
