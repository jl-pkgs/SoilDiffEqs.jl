using SoilDifferentialEquations, Test
using HydroTools
using OrdinaryDiffEqTsit5
using Plots
gr(framestyle=:box, legend=:topright)

includet("main_pkgs.jl")



using Dates
## 加入蒸发的信号
dates = DateTime(2010, 6):Hour(1):DateTime(2010, 9, 30)
ntime = length(dates)
lon, lat = 120.0, 30.0

doys = dayofyear.(dates)
hours = hour.(dates)

Rsi = map(t -> cal_Rsi_toa_hour(hours[t]; lat, J=doys[t]), eachindex(dates)) # MJ per period
Rn = Rsi / 4

# MJ h-1
Ta = 20.0
Δ = cal_slope(Ta)
γ = cal_gamma(Ta)
ET = Rn * Δ / (Δ + γ) * 1

## 从第一层土壤进行蒸发

## 将ET划分到土壤的每一层。
begin
  Δ = 0.1
  dt = 3600
  N = floor(Int, 2 / Δ)
  dz = fill(Δ, N) # 2m
  θ = LinRange(0.36, 0.2, N) |> collect

  soil = init_soil(; dt, zwt=-2.5)
  cal_θEψE!(soil) # update θE, ψE
  soil.θ = soil.θE[1:N]

  solver = Tsit5()
  @time R = solve_SM_Zeng2009(soil; solver, verbose=true, ET)
end


## 设置根系分布情况
"""
- z: in cm, 向下为正
"""
root_fraction(z_cm::AbstractVector; β) = β .^ z_cm # 

root_fraction(soil::Soil; β=0.943) = root_fraction(-soil.z_cm[0:N]; β)


begin
  βs = [0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97]

  p = plot(;legendposition=:bottomright, ylabel="Depths (m)", xlabel="Root Fraction (%)")
  for β in βs
    root = root_fraction(soil; β)
    plot!(p, -diff(root), soil.z[1:N], label="β =$β")
  end
  p
end






begin
  z = soil.z[1:N]
  p = plot(; size=(1000, 600))
  days = [1, 2, 3, 5, 7, 20, 30, 60, 90]
  for day in days
    t = day * 24
    plot!(p, R[t, :], z, label="day=$day")
  end
  plot!(p, soil.θE[1:N], z, label="θE", color=:black, linewidth=2)
  p
end


## 不同深度的时间序列
begin
  layers = [1, 2, 3, 4, 5, 7, 10]
  p = plot()
  t = 1:24*90
  for l in layers
    depth = round(Int, l * Δ * 100)
    plot!(p, dates[t], R[t, l], label="$(depth)cm")
  end
  p
end

## 绘制均衡水位
# zwt不变的情况下，均衡水位恒定
