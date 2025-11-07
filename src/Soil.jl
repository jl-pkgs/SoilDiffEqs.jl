export Soil


# @with_kw_noshow mutable struct SoilThermal{FT}
# end

# @with_kw_noshow mutable struct SoilHydraulic{FT,P<:AbstractSoilParam{FT}}
# end



@with_kw_noshow mutable struct Soil{FT,P<:AbstractSoilParam{FT}}
  N::Int = 10                        # layers of soil
  ibeg::Int = 1                      # index of the first layer，边界层条件指定
  inds_obs::Vector{Int} = ibeg:N     # indices of observed layers

  dt::Float64 = 3600.0               # 时间步长, seconds
  z::OffsetVector{FT} = zeros(FT, N + 2) # m, 向下为负
  Δz₊ₕ::Vector{FT} = zeros(FT, N + 1)
  z₊ₕ::Vector{FT} = zeros(FT, N + 1)
  Δz::Vector{FT} = zeros(FT, N + 1)

  z_cm::OffsetVector{FT} = z * 100      # cm, 向下为负
  Δz₊ₕ_cm::Vector{FT} = Δz₊ₕ * 100
  Δz_cm::Vector{FT} = Δz * 100

  # 水分
  θ::Vector{FT} = fill(0.1, N + 1)      # θ [m3 m-3]
  Q::Vector{FT} = zeros(FT, N)          # [cm h-1]
  K::Vector{FT} = zeros(FT, N + 1)      # hydraulic conductivity，[cm h-1], 1:N+1
  K₊ₕ::Vector{FT} = zeros(FT, N + 1)    # hydraulic conductivity at interface, [cm h-1]
  ∂θ∂ψ::Vector{FT} = zeros(FT, N + 1)   # specific moisture capacity, dθ/dΨ, [cm-1], 临时变量
  ψ::Vector{FT} = zeros(FT, N + 1)      # [cm]，约干越负
  ψ_next::Vector{FT} = zeros(FT, N + 1) # ψ[N+1/2], [cm], 临时变量
  θ0::FT = FT(0.0)                      # [m3 m-3]
  ψ0::FT = FT(0.0)                      # [cm]
  Q0::FT = FT(0.0)                      # [cm h-1] 下渗速率，向下为负
  sink::Vector{FT} = fill(0.0, N + 1)   # 蒸发项, [cm per unit time]
  θ_prev::Vector{FT} = zeros(FT, N + 1) # backup of θ  
  ψ_prev::Vector{FT} = zeros(FT, N + 1) # backup of ψ

  θE::Vector{FT} = zeros(FT, N + 1)     # equilibrium water content, [m3 m-3]
  ψE::Vector{FT} = zeros(FT, N + 1)     # water potential, [cm]

  # 地下水
  zwt::FT = FT(0.0)                  # groundwater depth, [m], 为了与z单位一致
  jwt::Int = 0                       # non-groundwater index
  wa::FT = FT(5000.0)                # water amount in aquifer, [mm]，潜水含水层
  uex::FT = FT(0.0)                  # 超出地表的水量, [mm], [kg m-2] 以地表径流的形式排放
  recharge::FT = FT(0.0)             # recharge rate, [cm/h]
  drainage::FT = FT(0.0)             # drainage rate, [cm/h]
  Sy_r::Vector{FT} = fill(0.2, N + 1)  # recharge, specific yield, [m3 m-3]
  Sy_d::Vector{FT} = fill(0.2, N + 1)  # discharge, specific yield, [m3 m-3]
  Sy_e::Vector{FT} = fill(0.2, N + 1)  # 均衡状态, specific yield, [m3 m-3]

  # 温度
  Tsoil::Vector{FT} = fill(NaN, N) # [°C]
  κ₊ₕ::Vector{FT} = zeros(FT, N)    # thermal conductivity at interface [W m-1 K-1]
  F::Vector{FT} = zeros(FT, N)     # heat flux, [W m-2]
  Tsurf::FT = FT(NaN)                  # surface temperature, [°C]
  F0::FT = FT(NaN)                     # heat flux at the surface, [W m-2]，向下为负
  G::FT = FT(NaN)                      # [W m-2]，土壤热通量

  ## 植被参数
  f_root::Vector{FT} = zeros(FT, N) # 根系分布比例
  βi::Vector{FT} = zeros(FT, N)     # 每层土壤水限制, dryness
  w::Vector{FT} = zeros(FT, N)      # 土壤水限制权重，f_root * βi
  β::FT = FT(0.1)                   # 总的土壤水限制

  ## Parameter: [水力] + [热力]参数
  method_retention::String = "van_Genuchten"
  param::SoilParam{FT,P} = SoilParam{FT,P}(; N, method_retention)

  # ODE求解临时变量
  u::Vector{FT} = fill(NaN, N + 1)  # [°C], 为了从ibeg求解地温，定义的临时变量
  du::Vector{FT} = fill(NaN, N + 1) # [°C]

  # 三角阵求解临时变量
  a::Vector{FT} = zeros(FT, N + 1)
  b::Vector{FT} = zeros(FT, N + 1)
  c::Vector{FT} = zeros(FT, N + 1)
  d::Vector{FT} = zeros(FT, N + 1)
  e::Vector{FT} = zeros(FT, N + 1)
  f::Vector{FT} = zeros(FT, N + 1)

  timestep::Int = 0                  # 迭代次数
end


function Soil{FT}(; method_retention::String="van_Genuchten", kw...) where {FT<:Real}
  if method_retention == "van_Genuchten"
    P = VanGenuchten{FT}
  elseif method_retention == "Campbell"
    P = Campbell{FT}
  end
  Soil{FT,P}(; method_retention, kw...)
end

function Soil(Δz::Vector{FT}; kw...) where {FT}
  N = length(Δz)
  Δz2 = zeros(FT, N+1)
  Δz2[1:N] = Δz

  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz=Δz2, Δz₊ₕ, kw...)
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


function Base.show(io::IO, x::Soil{T}) where {T<:Real}
  param = x.param

  printstyled(io, "Soil{$T}: ", color=:blue)
  printstyled(io, "N = $(x.N), ibeg=$(x.ibeg), ", color=:blue, underline=true)
  print_index(io, x.inds_obs; prefix="inds_obs =")

  printstyled(io, "Vegetation: \n", color=:blue, bold=true)
  print_var(io, x, :f_root)

  printstyled(io, "Soil Temperature: \n", color=:blue, bold=true)
  print_var(io, x, :Tsoil)
  print_var(io, x, :Tsurf)

  printstyled(io, "Soil Moisture: \n", color=:blue, bold=true)
  print_var(io, x, :K)
  print_var(io, x, :ψ)
  print_var(io, x, :θ)
  print_var(io, x, :sink)
  print_var(io, x, :Q)

  str = @sprintf(" - Q0 = %.3f, θ0 = %.3f, ψ0 = %.3f \n", x.Q0, x.θ0, x.ψ0)
  printstyled(io, str, color=:white, bold=false)

  # groundwater
  print_var(io, x, :zwt)
  print_var(io, x, :wa)

  # printstyled(io, "param_water: ", color=:blue, bold=true)
  # show(io, x.param_water)
  show(io, param)
  return nothing
end

using OffsetArrays

"""
    soil_depth_init(Δz::AbstractVector)
    
Soil depth initialization

```julia
z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
```
"""
function soil_depth_init(Δz::AbstractVector)
  # Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
  # z_{i+1/2}
  N = length(Δz)
  z = OffsetArray(zeros(N + 2), 0:N+1)    # 虚拟层N+1
  # Δz₊ₕ = OffsetArray(zeros(N + 1), 0:N)
  Δz₊ₕ = zeros(N + 1)
  z₊ₕ = zeros(N + 1)
  z₋ₕ = zeros(N + 1)

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
  # Δz₊ₕ[0] = 0.5 * Δz[1]
  for i = 1:N-1
    Δz₊ₕ[i] = z[i] - z[i+1]
  end
  Δz₊ₕ[N] = 0.5 * Δz[N]

  ## z₋ₕ
  z₋ₕ[1] = 0
  z₋ₕ[2:end] = z₊ₕ[1:end-1]
  (; z, z₋ₕ, z₊ₕ, Δz₊ₕ)
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
# z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
# cal_Δz(z₊ₕ) == Δz
# cal_Δz₊ₕ(z, z₊ₕ) == Δz₊ₕ
