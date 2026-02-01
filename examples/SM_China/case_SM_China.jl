using SoilDifferentialEquations, Ipaper, RTableTools
using Printf, Dates
import ModelParams: sceua
include("main_plot.jl")

# ========== Data Loading ==========
d = fread(joinpath(@__DIR__, "SM_J1193.csv"))
dates = d[:, 1]
_A = d[:, 2:end] ./ 100 |> Matrix |> drop_missing

function interp_depth(x, A, xout)
  ntime = size(A, 1)
  yout = zeros(ntime, length(xout))
  for i in 1:ntime
    yout[i, :] .= approx(x, view(A, i, :), xout)
  end
  yout
end

A = interp_depth([10, 20, 30, 40, 50, 60, 80, 100], _A, 10:10:100.)

# ========== Model Functions ==========
function model_sim(theta)
  soil = init_soil(; θ0, soil_type=8)
  SM_UpdateParam!(soil, theta)
  method = options.method_solve
  method == "Bonan" ? solve_SM_Bonan(soil, options.θ_surf) : 
    solve_SM_ODE(soil, options.θ_surf; solver=Tsit5())
end

function goal(theta)
  yobs = options.yobs
  ysim = model_sim(theta)
  ncol = size(yobs, 2)
  n = ncol - ibeg + 1
  sum(i -> -of_NSE(view(yobs, :, i), view(ysim, :, i)), ibeg:ncol) / n
end


# ========== Init Soil ==========
function init_soil(; θ0=0.3, dt=3600.0, soil_type=7)
  Δz = [5.0; repeat([10], 10)] ./ 100
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  (; method_retention, same_layer, ibeg) = options
  
  θ = fill(0.2, N)
  θ[2:end] .= θ0
  θ[1] = θ0[1]
  
  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par; use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
  cal_ψ!(soil)
  cal_K!(soil)
  soil
end

# ========== Main ==========
# 配置：10cm作为边界层，从20cm开始模拟
surface_depth_cm = 10
target_depths = collect(10:10:100.)
surface_idx = findfirst(==(surface_depth_cm), target_depths)
ibeg = surface_idx + 2  # 自动计算：从surface的下一层开始

θ_surf = A[:, surface_idx]
θ0 = A[1, :]
yobs = A[:, ibeg-1:end]

set_option!(; method_retention="van_Genuchten", yobs, θ_surf, ibeg, same_layer=false)

soil = init_soil(; θ0, soil_type=7)
lower, upper = SM_paramBound(soil)
theta0 = SM_param2theta(soil)

println("Initial NSE: $(-goal(theta0))")
ysim0 = model_sim(theta0)
plot_result(; ysim=ysim0, yobs, dates, depths=target_depths, ibeg, filename=joinpath(@__DIR__, "plot_initial.png"))

# Optimize
@time theta, feval, exitflag = sceua(goal, theta0, lower, upper; maxn=10_000)
println("Optimized. feval=$feval")

serialize(joinpath(@__DIR__, "theta"), theta)
theta = deserialize(joinpath(@__DIR__, "theta"))

SM_UpdateParam!(soil, theta)
ysim_opt = model_sim(theta)
plot_result(; ysim=ysim_opt, yobs, dates, depths=target_depths, ibeg, filename=joinpath(@__DIR__, "plot_optimized.png"))
