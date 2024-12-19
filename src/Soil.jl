export SoilParam, Soil
using Parameters
using Printf

export Soil, SoilParam, ParamVanGenuchten

## 结构体形式的参数
abstract type AbstractSoilParam{FT} end

@with_kw mutable struct ParamVanGenuchten{T} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  θ_res::T = 0.075       # [m3 m-3]
  Ksat::T = 34 / 3600    # [cm s-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0 - 1.0 / n
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
  m::Vector{FT} = fill(0.5, N)         # [-]，优化时的可选参数

  ψ_sat::Vector{FT} = fill(-10.0, N)   # [cm]
  b::Vector{FT} = fill(4.0, N)         # [-]

  ## Parameter: 土壤热力
  κ::Vector{FT} = fill(2.0, N)         # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = fill(2.0 * 1e6, N)  # volumetric heat capacity [J m-3 K-1]
end


@with_kw_noshow mutable struct Soil{FT}
  N::Int = 10                        # layers of soil
  ibeg::Int = 1                      # index of the first layer，边界层条件指定
  inds_obs::Vector{Int} = ibeg:N     # indices of observed layers

  dt::Float64 = 3600                 # 时间步长, seconds
  z::Vector{FT} = zeros(FT, N)       # m, 向下为负
  z₊ₕ::Vector{FT} = zeros(FT, N)
  Δz::Vector{FT} = zeros(FT, N)
  Δz₊ₕ::Vector{FT} = zeros(FT, N)

  z_cm::Vector{FT} = z * 100         # cm, 向下为负
  Δz_cm::Vector{FT} = Δz * 100
  Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100

  # 水分
  θ::Vector{FT} = fill(0.1, N)       # θ [m3 m-3]
  Q::Vector{FT} = zeros(FT, N)       # [cm/s]
  K::Vector{FT} = zeros(FT, N)       # hydraulic conductivity，[cm/s]
  K₊ₕ::Vector{FT} = zeros(FT, N - 1)  # hydraulic conductivity at interface, [cm/s]
  Cap::Vector{FT} = zeros(FT, N)     # specific moisture capacity, dθ/dΨ, [cm-1], 临时变量
  ψ::Vector{FT} = zeros(FT, N)       # [cm]，约干越负
  ψ_next::Vector{FT} = zeros(FT, N)  # ψ[N+1/2], [cm], 临时变量
  θ0::FT = FT(0.0)                   # [m3 m-3]
  ψ0::FT = FT(0.0)                   # [cm]
  Q0::FT = FT(0.0)                   # [cm/s] 下渗速率，向下为负
  sink::Vector{FT} = fill(0.0, N)    # 蒸发项, [cm per unit time]
  θ_prev::Vector{FT} = zeros(FT, N)  # backup of θ
  ψ_prev::Vector{FT} = zeros(FT, N)  # backup of ψ

  # 地下水
  zwt::FT = FT(0.0)                  # groundwater depth, [m], 为了与z单位一致
  wa::FT = FT(5000.0)                # water amount in aquifer, [mm]，潜水含水层
  uex::FT = FT(0.0)                  # 超出地表的水量, [mm], [kg m-2] 以地表径流的形式排放
  recharge::FT = FT(0.0)             # recharge rate, [mm/s]
  drainage::FT = FT(0.0)             # drainage rate, [mm/s]
  Sy::Vector{FT} = fill(0.02, N)     # specific yield, [m3 m-3]

  # 温度
  Tsoil::Vector{FT} = fill(NaN, N)   # [°C]
  κ₊ₕ::Vector{FT} = zeros(FT, N - 1)  # thermal conductivity at interface [W m-1 K-1]
  F::Vector{FT} = zeros(FT, N)       # heat flux, [W m-2]
  Tsurf::FT = FT(NaN)                # surface temperature, [°C]
  F0::FT = FT(NaN)                   # heat flux at the surface, [W m-2]，向下为负
  G::FT = FT(NaN)                    # [W m-2]，土壤热通量

  ## Parameter: [水力] + [热力]参数
  param::SoilParam{FT} = SoilParam{FT}(; N)
  param_water::ParamVanGenuchten{FT} = ParamVanGenuchten{FT}()

  # ODE求解临时变量
  u::Vector{FT} = fill(NaN, N)  # [°C], 为了从ibeg求解地温，定义的临时变量
  du::Vector{FT} = fill(NaN, N) # [°C]

  # 三角阵求解临时变量
  a::Vector{FT} = zeros(FT, N)
  b::Vector{FT} = zeros(FT, N)
  c::Vector{FT} = zeros(FT, N)
  d::Vector{FT} = zeros(FT, N)
  e::Vector{FT} = zeros(FT, N)
  f::Vector{FT} = zeros(FT, N)

  timestep::Int = 0                  # 迭代次数
end

function Soil(Δz::Vector{FT}; kw...) where {FT}
  N = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, kw...)
  # update K and ψ
  cal_K!(soil)
  cal_ψ!(soil)
  return soil
end

# θ = fill(0.1, N)
# ψ = van_Genuchten_ψ.(θ; param=param_water)
# θ0 = 0.267
# ψ0 = van_Genuchten_ψ(θ0; param=param_water)
# dt = 5 # [s]
# sink = ones(N) * 0.3 / 86400 # [cm s⁻¹], 蒸发速率

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


function Base.show(io::IO, x::Soil{T}) where {T<:Real}
  param = x.param

  printstyled(io, "Soil{$T}: ", color=:blue)
  printstyled(io, "N = $(x.N), ibeg=$(x.ibeg), ", color=:blue, underline=true)
  print_index(io, x.inds_obs; prefix="inds_obs =")

  printstyled(io, "Soil Temperature: \n", color=:blue, bold=true)
  print_var(io, x, :Tsoil)
  print_var(io, x, :Tsurf)

  printstyled(io, "Soil Moisture: \n", color=:blue, bold=true)
  print_var(io, x, :K, scale=1 / 3600) # [cm/s] to [cm h-1]
  print_var(io, x, :ψ)
  print_var(io, x, :θ)
  print_var(io, x, :sink)
  print_var(io, x, :θ0)
  print_var(io, x, :ψ0)

  # groundwater
  print_var(io, x, :zwt)
  print_var(io, x, :wa)

  # printstyled(io, "param_water: ", color=:blue, bold=true)
  # show(io, x.param_water)
  show(io, param)
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


"""
    soil_depth_init(Δz::AbstractVector)
    
Soil depth initialization

```julia
z, z₊ₕ, dz₊ₕ = soil_depth_init(Δz)
```
"""
function soil_depth_init(Δz::AbstractVector)
  # Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
  # z_{i+1/2}
  N = length(Δz)

  z = zeros(N)
  z₊ₕ = zeros(N)
  dz₊ₕ = zeros(N)

  z₊ₕ[1] = -Δz[1]
  for i = 2:N
    z₊ₕ[i] = z₊ₕ[i-1] - Δz[i] # on the edge
  end

  # Soil depth (m) at center of layer i (negative distance from surface)
  z[1] = 0.5 * z₊ₕ[1]
  for i = 2:N
    z[i] = 0.5 * (z₊ₕ[i-1] + z₊ₕ[i]) # on the center
  end

  # Thickness between between z(i) and z(i+1)
  for i = 1:N-1
    dz₊ₕ[i] = z[i] - z[i+1]
  end
  dz₊ₕ[N] = 0.5 * Δz[N]

  (; z, z₊ₕ, dz₊ₕ)
end


function cal_Δz(z)
  N = length(z)
  z₊ₕ = zeros(N)
  Δz = zeros(N)
  Δz[1] = 0 - z[1] * 2
  z₊ₕ[1] = -Δz[1]

  for i in 2:N
    Δz[i] = (z₊ₕ[i-1] - z[i]) * 2
    z₊ₕ[i] = z₊ₕ[i-1] - Δz[i]
  end
  Δz
end

# "face to center"
# function C2F(z::AbstractVector)
#   N = length(z)
#   z₊ₕ = zeros(N)
#   d = zeros(N)

#   z₊ₕ[1] = z[1] * 2
#   d[1] = z[1] * 2

#   @inbounds for i = 2:N
#     d[i] = 2(z[i] - z[i-1]) - d[i-1]
#     z₊ₕ[i] = z₊ₕ[i-1] + d[i]
#   end
#   z₊ₕ
# end

# "center to face"
# function F2C(z₊ₕ::AbstractVector)
#   N = length(z₊ₕ)
#   z = zeros(N)
#   z[1] = z₊ₕ[1] / 2
#   @inbounds for i = 2:N
#     z[i] = 0.5 * (z₊ₕ[i] + z₊ₕ[i-1])
#   end
#   z
# end
# export cal_Δz₊ₕ, C2F, F2C
export cal_Δz

# Δz = [2, 4, 6, 10]
# z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
# cal_Δz(z₊ₕ) == Δz
# cal_Δz₊ₕ(z, z₊ₕ) == Δz₊ₕ
