export AbstractSoilParam, ParamVanGenuchten, ParamCampbell
export SoilParam, get_soilpar

using Parameters
abstract type AbstractSoilParam{T<:Real} end

@with_kw mutable struct ParamVanGenuchten{T<:Real} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  θ_res::T = 0.075       # [m3 m-3]
  Ksat::T = 34.0         # [cm h-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0 - 1.0 / n
end

@with_kw mutable struct ParamCampbell{T<:Real} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  # θ_res::T = 0.075     # [m3 m-3]
  ψ_sat::T = -10.0       # [cm]
  Ksat::T = 34.0         # [cm h-1]
  b::T = 4.0             # [-]
end

function build_param(; method_retention::String="van_Genuchten", use_m::Bool=false,
  θ_sat::V, θ_res::V, Ksat::V, α::V, n::V, m::V, ψ_sat::V, b::V) where {V<:AbstractVector{<:Real}}
  T = eltype(θ_sat)
  if method_retention == "van_Genuchten"
    _m = use_m ? m : T(1) .- T(1) ./ n
    ψ_sat .= T(0) # update 20250517
    return ParamVanGenuchten{T}.(θ_sat, θ_res, Ksat, α, n, _m)
  elseif method_retention == "Campbell"
    return ParamCampbell{T}.(θ_sat, ψ_sat, Ksat, b)
  end
end

# 参数优化过程中，可能需要优化的参数
# 一个重要的经验教训，不要去优化`m`，NSE会下降0.2
@with_kw mutable struct SoilParam{FT, P<:AbstractSoilParam{FT}}
  ## Parameter: 土壤水力
  N::Int = 10
  method_retention::String = "van_Genuchten"     # "van_Genuchten" or "Campbell"
  use_m::Bool = false
  same_layer = false

  θ_sat::Vector{FT} = fill(FT(0.4), N)     # saturated water content, [m3 m-3]
  θ_res::Vector{FT} = fill(FT(0.1), N)     # residual water content, [m3 m-3]
  θ_fc::Vector{FT} = fill(FT(0.2), N)      # field capacity, [m3 m-3]
  Ksat::Vector{FT} = fill(FT(2.0), N)      # saturated hydraulic conductivity, [cm h-1]
  α::Vector{FT} = fill(FT(0.01), N)        # [cm-1]
  n::Vector{FT} = fill(FT(2.0), N)         # [-]
  m::Vector{FT} = fill(FT(0.5), N)         # [-]，优化时的可选参数，不建议优化

  ψ_sat::Vector{FT} = fill(FT(-10.0), N)   # [cm]
  b::Vector{FT} = fill(FT(4.0), N)         # [-]

  # soil moisture parameters
  param::Vector{P} = build_param(; method_retention, use_m, θ_sat, θ_res, Ksat, α, n, m, ψ_sat, b)
  ## Parameter: 土壤热力
  κ::Vector{FT} = fill(FT(2.0), N)         # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = fill(FT(2.0 * 1e6), N)  # volumetric heat capacity [J m-3 K-1]
end

function SoilParam{FT}(; method_retention::String="van_Genuchten", kw...) where {FT<:Real}
  if method_retention == "van_Genuchten"
    P = ParamVanGenuchten{FT}
    param = SoilParam{FT,P}(; method_retention, kw...)
    param.ψ_sat .= FT(0) # update 20250517
  elseif method_retention == "Campbell"
    P = ParamCampbell{FT}
    param = SoilParam{FT,P}(; method_retention, kw...)
  end
  return param
end


function SoilParam(N::Int, par::ParamCampbell{T};
  same_layer::Bool=true, kw...) where {T<:Real}
  
  (; θ_sat, ψ_sat, Ksat, b) = par
  SoilParam{T,ParamCampbell{T}}(; N,
    θ_sat=fill(θ_sat, N),
    ψ_sat=fill(ψ_sat, N),
    Ksat=fill(Ksat, N),
    b=fill(b, N),
    method_retention="Campbell", same_layer, kw...)
end

function SoilParam(N::Int, par::ParamVanGenuchten{T};
  use_m::Bool=false, 
  same_layer::Bool=true, kw...) where {T<:Real}

  (; θ_sat, θ_res, Ksat, α, n, m) = par
  _m = use_m ? m : (1 - 1 / n)
  SoilParam{T,ParamVanGenuchten{T}}(; N,
    ψ_sat = fill(T(0), N), # update 20250517
    θ_sat=fill(θ_sat, N),
    θ_res=fill(θ_res, N),
    Ksat=fill(Ksat, N),
    α=fill(α, N),
    n=fill(n, N),
    m=fill(_m, N),
    method_retention="van_Genuchten", same_layer, use_m, kw...)
end



function Update_SoilParam_Param!(soilparam::SoilParam{T}) where {T<:Real}
  (; N, method_retention, param, use_m) = soilparam
  (; θ_sat, θ_res, Ksat, α, m, n, ψ_sat, b) = soilparam

  if method_retention == "Campbell"
    for i in 1:N
      par = param[i]
      par.θ_sat = θ_sat[i]
      par.Ksat = Ksat[i]
      par.ψ_sat = ψ_sat[i]
      par.b = b[i]
    end
  elseif method_retention == "van_Genuchten"
    for i in 1:N
      par = param[i]
      par.θ_sat = θ_sat[i]
      par.θ_res = θ_res[i]
      par.Ksat = Ksat[i]
      par.α = α[i]
      par.n = n[i]
      par.m = use_m ? m[i] : (T(1) .- T(1) ./ n[i])
    end
  end
end

# 使用`StructVector`参数的好处：
# - 调用土壤水力函数时，可自动识别应调用的函数
function Base.show(io::IO, param::SoilParam{T}) where {T<:Real}
  (; use_m, same_layer) = param
  printstyled(io, "Parameters: \n", color=:blue, bold=true)
  println(typeof(param))
  # println(P)
  # println("[use_m = $use_m, same_layer = $same_layer]")

  println(io, "-----------------------------")
  print_var(io, param, :κ)
  print_var(io, param, :cv; scale=1e6)
  println(io, "-----------------------------")

  method_retention = param.method_retention
  subfix = same_layer ? " * 1" : " * N"
  np = use_m ? 6 : 5
  print_var(io, param, :θ_fc)

  print_selected(io, "van_Genuchten ($(np)p$subfix)", method_retention)
  print_var(io, param, :θ_sat)
  print_var(io, param, :θ_res)
  print_var(io, param, :Ksat)
  print_var(io, param, :α)
  print_var(io, param, :n)
  use_m && print_var(io, param, :m; used=use_m)
  print_selected(io, "Campbell (4p$subfix)", method_retention)
  printstyled(io, " - θ_sat, Ksat \n", color=:blue)

  print_var(io, param, :ψ_sat)
  print_var(io, param, :b)
  return nothing
end


function print_selected(io::IO, name::String, method_retention::String)
  if name[1:5] == method_retention[1:5]
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
