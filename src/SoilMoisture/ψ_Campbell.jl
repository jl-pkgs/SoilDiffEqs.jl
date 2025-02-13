"""
    Campbell(ψ::T, par::ParamCampbell{T})

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
@inline function Campbell(ψ::T, par::ParamCampbell{T}) where {T<:Real}
  θ = Campbell_θ(ψ, par)
  K = Campbell_K(θ, par)
  ∂θ∂ψ = Campbell_∂θ∂ψ(ψ, par)
  θ, K, ∂θ∂ψ
end

"""
    Campbell_θ(ψ, ψ_sat, θ_sat, b)
"""
@inline @fastmath function Campbell_θ(ψ::T, par::ParamCampbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, b) = par
  ψ <= ψ_sat ? θ_sat * (ψ / ψ_sat)^(-1 / b) : θ_sat
end

@inline @fastmath function Campbell_∂θ∂ψ(ψ::T, par::ParamCampbell{T}) where {T<:Real}
  (; ψ_sat, θ_sat, b) = par
  ψ <= ψ_sat ? -θ_sat / (b * ψ_sat) * (ψ / ψ_sat)^(-1 / b - 1) : 0.0
end

"""
    Campbell_K(θ, θ_sat, Ksat, b)
"""
@inline @fastmath function Campbell_K(θ::T, par::ParamCampbell{T}) where {T<:Real}
  (; θ_sat, Ksat, b) = par
  se = clamp(θ / θ_sat, 0.0, 1.0)
  Ksat * se^(2b + 3)
end

# Campbell 1974, Bonan 2019 Table 8.2
"""
    Campbell_ψ(θ, θ_sat, ψ_sat, b)
"""
@inline @fastmath function Campbell_ψ(θ::T, par::ParamCampbell{T}) where {T<:Real}
  (; θ_sat, ψ_sat, b) = par
  se = clamp(θ / θ_sat, 0.0, 1.0)
  ψ = ψ_sat * se^(-b)
  min(ψ, ψ_sat) # ψ为负值
end


Retention(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell(ψ, par)
Retention_K(θ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_K(θ, par)
Retention_θ(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_θ(ψ, par)
Retention_ψ(θ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_ψ(θ, par)
Retention_∂θ∂ψ(ψ::T, par::ParamCampbell{T}) where {T<:Real} = Campbell_∂θ∂ψ(ψ, par)

Retention(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention(ψ, par)
Retention_K(θ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_K(θ, par)
Retention_θ(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_θ(ψ, par)
Retention_ψ(θ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_ψ(θ, par)
Retention_∂θ∂ψ(ψ::T; par::AbstractSoilParam{T}) where {T<:Real} = Retention_∂θ∂ψ(ψ, par)

export Campbell, Campbell_ψ, Campbell_θ, Campbell_K
export Retention, Retention_K, Retention_θ, Retention_ψ, Retention_∂θ∂ψ
