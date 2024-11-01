export SoilParam, Soil
using Parameters
using Printf


# 参数优化过程中，可能需要优化的参数
# 一个重要的经验教训，不要去优化`m`，NSE会下降0.2
@with_kw mutable struct SoilParam{FT}
  ## Parameter: 土壤水力
  N::Int = 10
  method::String = "van_Genuchten"
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
  κ::Vector{FT} = fill(2.0, N)       # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = fill(2.0*1e6, N)      # volumetric heat capacity [J m-3 K-1]
end


@with_kw_noshow mutable struct Soil{FT}
  N::Int = 10                        # layers of soil
  ibeg::Int = 1                      # index of the first layer，边界层条件指定
  inds_obs::Vector{Int} = ibeg:N     # indices of observed layers

  dt::Float64 = 3600                 # 时间步长, seconds
  z::Vector{FT} = zeros(FT, N)       # cm, 向下为负
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

  ## Parameter: 土壤水力
  param_water::ParamVanGenuchten{FT} = ParamVanGenuchten{FT}()
  param::SoilParam{FT} = SoilParam{FT}(; N)

  # 温度
  Tsoil::Vector{FT} = fill(NaN, N)   # [°C]
  κ₊ₕ::Vector{FT} = zeros(FT, N - 1)  # thermal conductivity at interface [W m-1 K-1]
  F::Vector{FT} = zeros(FT, N)       # heat flux, [W m-2]
  TS0::FT = FT(NaN)                  # surface temperature, [°C]
  F0::FT = FT(NaN)                   # heat flux at the surface, [W m-2]，向下为负
  G::FT = FT(NaN)                    # [W m-2]，土壤热通量

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
  print_selected(io, "Cambell (4p$subfix)", method)
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
  print_var(io, x, :TS0)

  printstyled(io, "Soil Moisture: \n", color=:blue, bold=true)
  print_var(io, x, :K)
  print_var(io, x, :ψ)
  print_var(io, x, :θ)
  print_var(io, x, :sink)
  print_var(io, x, :θ0)
  print_var(io, x, :ψ0)

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
