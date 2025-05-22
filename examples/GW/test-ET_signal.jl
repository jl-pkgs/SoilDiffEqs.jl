using SoilDifferentialEquations, Test
using HydroTools
using OrdinaryDiffEqTsit5
using Plots, Dates
# gr(framestyle=:box, legend=:bottomleft)
pyplot(framestyle=:box, legend=:bottomleft)

include("main_pkgs.jl")
include("main_vis.jl")

function Init_VegParam!(soil; β=0.94)
  soil.f_root = root_fraction(soil; β)
end

# ylabel = "Depth (m)"
# plot(
#   plot(soil.f_root, soil.z[1:N]; label="", ylabel, xlabel="Root Fraction"), 
#   plot(cumsum(soil.f_root), soil.z[1:N]; label="", ylabel, xlabel="Accu Root Fraction"), 
#   size=(900, 500)
# )
# gr(framestyle=:box, )
# gr(display_type=:inline)

## 1. 将ET划分到土壤的每一层；模块已成。
begin
  # Δ = 0.1
  zwt = -2.5
  dt = 3600
  # N = floor(Int, 2 / Δ)
  # rep()
  dz = [5.0; repeat([10], 10); 50; 50] ./ 100 # 厚度
  # dz = fill(Δ, N) # 2m
  N = length(dz)
  θ = LinRange(0.36, 0.2, N) |> collect

  soil = init_soil(; dt, zwt)
  Init_VegParam!(soil)
  cal_θEψE!(soil) # update θE, ψE
  # soil.θ = soil.θE[1:N]
  soil.θ .= soil.param.θ_sat[1] .* 0.8

  W_prev = sum(soil.Δz[1:N] .* soil.θ[1:N]) * 1000 # [m] -> [mm]
  solver = Tsit5()

  ## 加入蒸发的信号
  dates = DateTime(2010, 4):Hour(1):DateTime(2010, 10, 31)
  ET = generate_ET(dates)
  @show sum(ET)

  # soil.zwt = -1.5
  soil.zwt = -2.5
  # soil.zwt = -5.0
  zwt = soil.zwt
  GW = false
  @time SM, SINK, Q, Z = solve_SM_Zeng2009(soil; solver, verbose=true, ET, reltol=1e-4, abstol=1e-4, GW)
  W_now = sum(soil.Δz[1:N] .* soil.θ[1:N]) * 1000 # [m] -> [mm]

  p1 = plot_soil_profile(soil)
  p2 = plot_depth_timeseries(soil)
  p3 = plot(dates, Z, title="(c) Groundwater Level (m)", label="", ylabel="Depth (m)")
  p4 = plot(dates, Q[:, end], title="(d) Discharge Rate [cm h-1]", label="")

  p = plot(p1, p2, p3, p4; size=(1600, 800), layout=(2, 2))
  savefig(p, "Figure1_ET_zwt=$(-zwt)_(GW=$GW).png")
  display(p)

  R = sum(Q[:, end] * 10) # 最后一层流向地下水 cm h-1, 水量不平衡
  (W_prev + R - W_now - sum(SINK) * 10) / W_prev * 100  # 误差 0.0125%
end

## 判断总的水量是否守恒

## m^3
## 严重干旱期间，出现了补给量大于排泄量的情况
begin
  k = 2
  q_in = Q[:, k-1]
  q_out = Q[:, k]
  et = SINK[:, k]

  p = plot()
  plot!(dates, -q_in, label="q_in")
  plot!(dates, q_out, label="q_out") # 正为补给，负为排泄
  dq = q_out - q_in
  plot!(dates, dq, label="dq")
  # plot!(dates, et, label="ET")
  p
end

begin
  k = 3
  q_in = Q[:, k-1]
  q_out = Q[:, k]
  et = SINK[:, k]

  p = plot()
  plot!(dates, -q_in, label="q_in")
  plot!(dates, q_out, label="q_out") # 正为补给，负为排泄
  dq = q_out - q_in
  plot!(dates, dq, label="dq")
  # plot!(dates, et, label="ET")
  p
end



## 绘制均衡水位
# zwt不变的情况下，均衡水位恒定
begin
  βs = [0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97]

  p = plot(; legendposition=:bottomleft,
    xlim=(0.6, 0.99), ylim=(-1, 0),
    ylabel="Depths (m)", xlabel="Accumulated Root Fraction (%)")
  vline!(p, [0.95], label=""; linestyle=:dash, color=:red)
  for β in βs
    root = root_fraction(soil; β)
    fraction = -diff(root)
    y = cumsum(fraction)
    plot!(p, y, soil.z[1:N], label="β =$β")
  end
  p
end
