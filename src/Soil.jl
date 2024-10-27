using Parameters
using Printf

# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = x^y
# @fastmath pow(x::Real, y::Real) = exp(y * log(x))

# Ksat: [cm/s]
abstract type AbstractSoilParam{FT} end

@with_kw mutable struct ParamVanGenuchten{T} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  θ_res::T = 0.075       # [m3 m-3]
  Ksat::T = 34 / 3600 # [cm s-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0
end

function Base.Vector(x::ParamVanGenuchten)
  (;θ_sat, θ_res, Ksat, α, n) = x
  [θ_sat, θ_res, Ksat, α, n]
end

# Bonan 2019, Table 8.3
function get_soilpar(soil_type::Int=1)
  soilparam = [
    # θ_sat, θ_res, α (cm⁻¹), n, Ksat (cm h⁻¹)
    0.38 0.068 0.008 1.09 0.2;   #  1,  Clay
    0.36 0.07  0.005 1.09 0.02;  #  2,  Silty  clay
    0.38 0.1   0.027 1.23 0.12;  #  3,  Sandy  clay
    0.41 0.095 0.019 1.31 0.26;  #  4,  Clay   loam
    0.43 0.089 0.01  1.23 0.07;  #  5,  Silty  clay loam
    0.39 0.1   0.059 1.48 1.31;  #  6,  Sandy  clay loam
    0.43 0.078 0.036 1.56 1.04;  #  7,  Loam
    0.45 0.067 0.02  1.41 0.45;  #  8,  Silty  loam
    0.41 0.065 0.075 1.89 4.42;  #  9,  Sandy  loam
    0.41 0.065 0.075 1.89 4.42;  #  10, Silty, no   data in Bonan2019
    0.41 0.057 0.124 2.28 14.59; #  11, Loamy  sand
    0.43 0.045 0.145 2.68 29.7   #  12, Sand
  ]
  θ_sat, θ_res, α, n, Ksat = soilparam[soil_type, :]
  Ksat = Ksat / 3600 # [cm h-1] to [cm s-1]
  ParamVanGenuchten(; θ_sat, θ_res, α, n, Ksat)
end

function get_soilpar(theta::AbstractVector)
  θ_sat, θ_res, Ksat, α, n = theta[1:5]
  ParamVanGenuchten(; θ_sat, θ_res, α, n, Ksat)
end


@with_kw_noshow mutable struct Soil{FT}
  n::Int = 10                        # layers of soil
  ibeg::Int = 1                      # index of the first layer，边界层条件指定
  inds_obs::Vector{Int} = ibeg:n     # indices of observed layers

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


function Base.show(io::IO, x::Soil{T}) where {T<:Real}
  printstyled(io, "Soil{$T}: ", color=:blue)
  printstyled(io, "n = $(x.n), ibeg=$(x.ibeg), ", color=:blue, underline=true)
  print_index(io, x.inds_obs; prefix = "inds_obs =")

  printstyled(io, "Soil Temperature: \n", color=:blue, bold=true)
  print_var(io, x, :κ)
  print_var(io, x, :cv; scale=1e6)
  print_var(io, x, :Tsoil)
  print_var(io, x, :TS0)

  printstyled(io, "Soil Moisture: \n", color=:blue, bold=true)
  print_var(io, x, :K)
  print_var(io, x, :ψ)
  print_var(io, x, :θ)
  print_var(io, x, :sink)
  print_var(io, x, :θ0)
  print_var(io, x, :ψ0)

  printstyled(io, "param_water: ", color=:blue, bold=true)
  show(io, x.param_water)
  return nothing
end

function print_var(io::IO, x, var; scale=nothing, digits=2)
  value = getfield(x, var)
  name = @sprintf("%-5s", string(var))
  printstyled(io, " - $name: ", color=:blue)
  if isnothing(scale)
    println(io, round.(value; digits))
  else
    println(io, "$(round.(x.cv/scale; digits)) * $scale")
  end
end

function print_index(io::IO, inds; prefix="", color=:blue, underline=true)
  if length(unique(diff(inds))) == 1
    n = length(inds)
    printstyled(io, "$prefix $(inds[1]):$(inds[end]) [n=$n] \n"; color, underline)
  else
    printstyled(io, "$prefix $inds \n"; color, underline)
  end
end

export get_soilpar
