"""
    Retention(ψ::T, par::ParamCampbell{T})

Campbell (1974) Retention relationships

# Usage
```julia
θ, K, ∂θ∂ψ = Retention(ψ, par)
θ = Retention_θ(ψ, par)
K = Retention_K(θ, par)
ψ = Retention_ψ(θ, par)
∂θ∂ψ = Retention_∂θ∂ψ(ψ, par)
```

# Arguments
+ `ψ`: Matric potential, cm
+ `param`
  - `ψ_sat`: Matric potential at saturation, [cm]
  - `θ_sat`: Volumetric water content at saturation
  - `b`    : Exponent
  - `Ksat`: Hydraulic conductivity at saturation [cm h-1]

# Examples
```julia
θ_sat = 0.25
ψ_sat = -25.0
b = 0.2
Ksat = 3.4e-03
θ, K, ∂θ∂ψ = Campbell(ψ, ψ_sat, θ_sat, Ksat, b)
```
"""
@fastmath function Retention(ψ::T, par::ParamCampbell{T}) where {T<:Real}
  θ = Retention_θ(ψ, par)
  K = Retention_K(θ, par)
  ∂θ∂ψ = Retention_∂θ∂ψ(ψ, par)
  θ, K, ∂θ∂ψ
end

"""
    Retention_θ(ψ, ψ_sat, θ_sat, b)
"""
@inline @fastmath function Retention_θ(ψ::T, par::ParamCampbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, b) = par
  ψ <= ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : θ_sat
end

@inline @fastmath function Retention_∂θ∂ψ(ψ::T, par::ParamCampbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, b) = par
  ψ <= ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : 0.0
end

"""
    Retention_K(θ, θ_sat, Ksat, b)
"""
@inline @fastmath function Retention_K(θ::T, par::ParamCampbell{T}) where {T<:Real}
  (; θ_sat, Ksat, b) = par
  se = clamp(θ / θ_sat, 0.0, 1.0)
  Ksat * se^(2b + 3)
end

# Campbell 1974, Bonan 2019 Table 8.2
"""
    Retention_ψ(θ, θ_sat, ψ_sat, b)
"""
@inline @fastmath function Retention_ψ(θ::T, par::ParamCampbell{T}) where {T<:Real}
  (; θ_sat, ψ_sat, b) = par
  se = clamp(θ / θ_sat, 0.0, 1.0)
  ψ = ψ_sat * se^(-b)
  min(ψ, ψ_sat) # ψ为负值
end


export Retention, Retention_θ, Retention_K, Retention_∂θ∂ψ, Retention_ψ
