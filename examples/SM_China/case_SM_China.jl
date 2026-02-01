using SoilDifferentialEquations
include("main_optim.jl")
include("main_plot.jl")
include("load_data.jl")

# dz = [2.5, 5, 5, 15, 45, 55]
# z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
# Δz = cal_Δz(z)
function init_soil(; θ0=0.3, dt=3600.0, soil_type=7)
  Δz = [5.0; repeat([10], 10)] ./ 100 # 厚度
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  (; method_retention, same_layer, ibeg) = options

  θ = fill(0.2, N)
  θ[2:end] .= θ0 # note here
  θ[1] = θ0[1]

  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par;
    use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
  cal_ψ!(soil)     # ψ, 以θ为导向
  cal_K!(soil)     # K
  soil
end

## 从第3层开始模拟
ibeg = 4                # 含有1个虚假的土壤层（5cm）, 第2层(10cm)是输入，第3层模拟(20cm)
θ_surf = A[:, ibeg-2]   # 表层SM时间序列
θ0 = A[1, :]            # 初始时刻的SM, 10cm~100cm
yobs = A[:, ibeg-1:end] # 观测值, 20cm~100cm
# inds_obs = [2:7; 9; 11] # 观测SM对应的层
set_option!(; method_retention="van_Genuchten", yobs, θ_surf, ibeg, same_layer=false)
options

soil = init_soil(; θ0, soil_type=7)
lower, upper = SM_paramBound(soil)
theta0 = SM_param2theta(soil)
ysim = model_sim(theta0;)

# ysim = model_sim(theta;) # theta undefined; run after optimization

goal(theta0)
plot_result(theta0)


begin
  f(theta) = goal(theta)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=10_000)
end

f_theta = "$(@__DIR__)/theta"
serialize(f_theta, theta)
theta = deserialize(f_theta)

SM_UpdateParam!(soil, theta)
plot_result(theta)
# Iteration =  11, nEvals = 10166, Best Cost = -0.91137
# 380.999683 seconds (10.90 M allocations: 48.204 GiB, 2.44% gc time, 0.41% compilation time: 100% of which was recompilation)
