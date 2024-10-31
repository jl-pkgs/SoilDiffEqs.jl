using Parameters
using Printf


@with_kw_noshow mutable struct Soil{FT}
  n::Int = 10                        # layers of soil
  ibeg::Int = 1                      # index of the first layer，边界层条件指定
  inds_obs::Vector{Int} = ibeg:n     # indices of observed layers

  dt::Float64 = 3600                 # 时间步长, seconds
  z::Vector{FT} = zeros(FT, n)       # cm, 向下为负
  z₊ₕ::Vector{FT} = zeros(FT, n)
  Δz::Vector{FT} = zeros(FT, n)
  Δz₊ₕ::Vector{FT} = zeros(FT, n)

  z_cm::Vector{FT} = z * 100         # cm, 向下为负
  Δz_cm::Vector{FT} = Δz * 100
  Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100

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

  ## Parameter: 土壤水力
  θ_sat::Vector{FT} = fill(0.4, n)     # saturated water content, [m3 m-3]
  θ_res::Vector{FT} = fill(0.1, n)     # residual water content, [m3 m-3]
  Ksat::Vector{FT} = fill(2.0 / 3600, n) # saturated hydraulic conductivity, [cm s-1], 2.0 cm/h
  α::Vector{FT} = fill(0.01, n)        # [m-1]
  n::Vector{FT} = fill(2.0, n)         # [-]
  m::Vector{FT} = fill(0.5, n)         # [-]，优化时的可选参数

  ψ_sat::Vector{FT} = fill(-10.0, n)   # [cm]
  b::Vector{FT} = fill(4.0, n)         # [-]

  # 温度
  Tsoil::Vector{FT} = fill(NaN, n)   # [°C]
  κ₊ₕ::Vector{FT} = zeros(FT, n - 1)  # thermal conductivity at interface [W m-1 K-1]
  F::Vector{FT} = zeros(FT, n)       # heat flux, [W m-2]
  TS0::FT = FT(NaN)                  # surface temperature, [°C]
  F0::FT = FT(NaN)                   # heat flux at the surface, [W m-2]，向下为负
  G::FT = FT(NaN)                    # [W m-2]，土壤热通量

  ## Parameter: 土壤热力
  κ::Vector{FT} = zeros(FT, n)       # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = zeros(FT, n)      # volumetric heat capacity [J m-3 K-1]

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
  print_index(io, x.inds_obs; prefix="inds_obs =")

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
