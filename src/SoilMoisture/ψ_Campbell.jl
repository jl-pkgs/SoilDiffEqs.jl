"""
    Campbell(ψ, ψ_sat, θ_sat, Ksat, b)

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
θ_sat = 0.25
ψ_sat = -25.0
b = 0.2
Ksat = 3.4e-03
θ, K, ∂θ∂ψ = Campbell(ψ, ψ_sat, θ_sat, Ksat, b)
```
"""
@fastmath function Campbell(ψ::T, ψ_sat::T, θ_sat::T, Ksat::T, b::T) where {T<:Real}
  θ = Campbell_θ(ψ, ψ_sat, θ_sat, b)
  K = Campbell_K(θ, θ_sat, Ksat, b)
  ∂θ∂ψ = Campbell_∂θ∂ψ(ψ, ψ_sat, θ_sat, b)
  θ, K, ∂θ∂ψ
end

"""
    Campbell_θ(ψ, ψ_sat, θ_sat, b)
"""
@inline @fastmath function Campbell_θ(ψ::T, ψ_sat::T, θ_sat::T, b::T) where {T<:Real}
  ψ <= ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : θ_sat
end

@inline @fastmath function Campbell_∂θ∂ψ(ψ::T, ψ_sat::T, θ_sat::T, b::T) where {T<:Real}
  ψ <= ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : 0.0
end

"""
    Campbell_K(θ, θ_sat, Ksat, b)
"""
@inline @fastmath function Campbell_K(θ::T, θ_sat::T, Ksat::T, b::T) where {T<:Real}
  se = clamp(θ / θ_sat, 0.0, 1.0)
  Ksat * se^(2b + 3)
end

# Campbell 1974, Bonan 2019 Table 8.2
"""
    Campbell_ψ(θ, θ_sat, ψ_sat, b)
"""
@inline @fastmath function Campbell_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  se = clamp(θ / θ_sat, 0.0, 1.0)
  ψ = ψ_sat * se^(-b)
  min(ψ, ψ_sat) # ψ为负值
end


export Campbell, Campbell_ψ, Campbell_θ, Campbell_K
