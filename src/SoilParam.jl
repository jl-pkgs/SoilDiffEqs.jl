export AbstractSoilParam, ParamVanGenuchten, ParamCampbell
export SoilParam, get_soilpar

abstract type AbstractSoilParam{FT} end

@with_kw mutable struct ParamVanGenuchten{T} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  θ_res::T = 0.075       # [m3 m-3]
  Ksat::T = 34 / 3600    # [cm s-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0 - 1.0 / n
end

@with_kw mutable struct ParamCampbell{T} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  # θ_res::T = 0.075     # [m3 m-3]
  ψ_sat::T = -10.0       # [cm]
  Ksat::T = 34 / 3600    # [cm s-1]
  b::T = 4.0             # [-]
end


function Base.Vector(x::ParamVanGenuchten)
  (; θ_sat, θ_res, Ksat, α, n, m) = x
  [θ_sat, θ_res, Ksat, α, n, m]
end
function Base.Vector(x::ParamCampbell)
  (; θ_sat, ψ_sat, Ksat, b) = x
  [θ_sat, ψ_sat, Ksat, b]
end
Base.collect(x::AbstractSoilParam) = Vector(x)


function get_soilpar(::Type{ParamVanGenuchten}, soil_type::Int=1)
  # Bonan 2019, Table 8.3
  soilparam = [
    # θ_sat, θ_res, α (cm⁻¹), n, Ksat (cm h⁻¹)
    0.38 0.068 0.008 1.09 0.2;   #  1, 11, Clay
    0.36 0.070 0.005 1.09 0.02;  #  2, 10, Silty clay
    0.38 0.100 0.027 1.23 0.12;  #  3,  9, Sandy clay
    0.41 0.095 0.019 1.31 0.26;  #  4,  8, Clay  loa
    0.43 0.089 0.010 1.23 0.07;  #  5,  7, Silty clay loam
    0.39 0.100 0.059 1.48 1.31;  #  6,  6, Sandy clay loam
    0.43 0.078 0.036 1.56 1.04;  #  7,  5, Loam
    0.45 0.067 0.020 1.41 0.45;  #  8,  4, Silty loam
    0.41 0.065 0.075 1.89 4.42;  #  9,  3, Sandy loam
    0.45 0.067 0.020 1.41 0.45;  #  10, NA, Silty, nodata, used Silty loam
    0.41 0.057 0.124 2.28 14.59; #  11, 2, Loamy sand
    0.43 0.045 0.145 2.68 29.7   #  12, 1, Sand
  ]
  θ_sat, θ_res, α, n, Ksat = soilparam[soil_type, :]
  Ksat = Ksat / 3600 # [cm h-1] to [cm s-1]
  ParamVanGenuchten(; θ_sat, θ_res, α, n, Ksat)
end

function get_soilpar(::Type{ParamVanGenuchten}, theta::AbstractVector)
  θ_sat, θ_res, Ksat, α, n = theta[1:5]
  ParamVanGenuchten(; θ_sat, θ_res, α, n, Ksat)
end

function get_soilpar(::Type{ParamCampbell}, soil_type::Int=1)
  # Bonan 2019, Table 8.3
  soilparam = [
    # θ_sat, ψ_sat (cm), b, Ksat (cm h⁻¹)
    0.482 -40.5 11.4 0.46;  # 11, Clay
    0.492 -49.0 10.4 0.37;  # 10, Silty clay
    0.426 -15.3 10.4 0.78;  # 9, Sandy clay
    0.476 -63.0 8.52 0.88;  # 8, Clay loam
    0.477 -35.6 7.75 0.61;  # 7, Silty clay loam
    0.420 -29.9 7.12 2.27;  # 6, Sandy clay loam
    0.451 -47.8 5.39 2.50;  # 5, Loam
    0.485 -78.6 5.30 2.59;  # 4, Silty loam
    0.435 -21.8 4.90 12.48; # 3, Sandy loam
    0.485 -78.6 5.30 2.59;  # NA, Silty, nodata, used Silty loam (4)
    0.410 -9.0 4.38 56.28;  # 2, Loamy sand
    0.395 -12.1 4.05 63.36  # 1, Sand
  ]
  θ_sat, ψ_sat, b, Ksat = soilparam[soil_type, :]
  Ksat = Ksat / 3600 # [cm h-1] to [cm s-1]
  ParamCampbell(; θ_sat, ψ_sat, b, Ksat)
end

function get_soilpar(::Type{ParamCampbell}, theta::AbstractVector)
  θ_sat, ψ_sat, b, Ksat = theta[1:4]
  ParamCampbell(; θ_sat, ψ_sat, b, Ksat)
end


const TypeRetention = Union{Type{ParamVanGenuchten},Type{ParamCampbell}}

function get_soilpar(soil_type::Int=1; type::TypeRetention=ParamVanGenuchten)
  get_soilpar(type, soil_type)
end

function get_soilpar(theta::AbstractVector; type::TypeRetention=ParamVanGenuchten)
  get_soilpar(type, theta)
end

# 参数优化过程中，可能需要优化的参数
# 一个重要的经验教训，不要去优化`m`，NSE会下降0.2
@with_kw mutable struct SoilParam{FT}
  ## Parameter: 土壤水力
  N::Int = 10
  method::String = "van_Genuchten"     # "van_Genuchten" or "Campbell"
  use_m::Bool = false
  same_layer = false
  
  θ_sat::Vector{FT} = fill(0.4, N)     # saturated water content, [m3 m-3]
  θ_res::Vector{FT} = fill(0.1, N)     # residual water content, [m3 m-3]
  Ksat::Vector{FT} = fill(2.0 / 3600, N) # saturated hydraulic conductivity, [cm s-1]
  α::Vector{FT} = fill(0.01, N)        # [m-1]
  n::Vector{FT} = fill(2.0, N)         # [-]
  m::Vector{FT} = fill(0.5, N)         # [-]，优化时的可选参数，不建议优化

  ψ_sat::Vector{FT} = fill(-10.0, N)   # [cm]
  b::Vector{FT} = fill(4.0, N)         # [-]
  
  # soil moisture parameters
  param::StructVector{<:AbstractSoilParam{FT}} = build_param(; method, θ_sat, θ_res, Ksat, α, n, m, ψ_sat, b)
  ## Parameter: 土壤热力
  κ::Vector{FT} = fill(2.0, N)         # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = fill(2.0 * 1e6, N)  # volumetric heat capacity [J m-3 K-1]
end

function build_param(; method::String="van_Genuchten",
  θ_sat::V, θ_res::V, Ksat::V, α::V, n::V, m::V, ψ_sat::V, b::V) where {V<:AbstractVector{<:Real}}
  FT = eltype(θ_sat)
  if method == "Campbell"
    return StructArray{ParamCampbell{FT}}(; θ_sat, ψ_sat, Ksat, b)
  elseif method == "van_Genuchten"
    return StructArray{ParamVanGenuchten{FT}}(; θ_sat, θ_res, Ksat, α, n, m)
  end
end

# 使用`StructVector`参数的好处：
# - 调用土壤水力函数时，可自动识别应调用的函数
function Base.show(io::IO, param::SoilParam{T}) where {T<:Real}
  (; use_m, same_layer) = param
  printstyled(io, "Parameters: \n", color=:blue, bold=true)
  # println("[use_m = $use_m, same_layer = $same_layer]")

  println(io, "-----------------------------")
  print_var(io, param, :κ)
  print_var(io, param, :cv; scale=1e6)
  println(io, "-----------------------------")

  method = param.method
  subfix = same_layer ? " * 1" : " * N"
  np = use_m ? 6 : 5
  print_selected(io, "van_Genuchten ($(np)p$subfix)", method)
  print_var(io, param, :θ_sat)
  print_var(io, param, :θ_res)
  print_var(io, param, :Ksat; scale=1e-3)
  print_var(io, param, :α)
  print_var(io, param, :n)
  use_m && print_var(io, param, :m; used=use_m)
  print_selected(io, "Campbell (4p$subfix)", method)
  printstyled(io, " - θ_sat, Ksat \n", color=:blue)

  print_var(io, param, :ψ_sat)
  print_var(io, param, :b)
  return nothing
end


function print_selected(io::IO, name::String, method::String)
  if name[1:5] == method[1:5]
    printstyled(io, "   [$name]\n", bold=true, color=:green)
  else
    printstyled(io, "   [$name]\n", bold=true)
  end
end

function print_var(io::IO, x, var; scale=nothing, digits=3, color=:blue, used=true)
  value = getfield(x, var)
  name = @sprintf("%-5s", string(var))
  _color = used ? color : :white
  printstyled(io, " - $name: "; color=_color)
  if isnothing(scale)
    println(io, round.(value; digits))
  else
    println(io, "$(round.(value/scale; digits)) * $scale")
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
