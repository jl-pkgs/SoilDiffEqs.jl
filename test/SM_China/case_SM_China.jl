using SoilDifferentialEquations
using SoilDifferentialEquations.GlobalOptions

GlobalOptions.options = Options()
options = GlobalOptions.options
options.method_retention = "van_Genuchten"
options


include("main_optim.jl")
include("main_plot.jl")
include("load_data.jl")

## 从第二层开始模拟
θ_surf = A[:, 1]   # 表层SM时间序列
θ0 = A[1, :]       # 初始时刻的SM
ibeg = 2           # 含有1个虚假的土壤层（5cm）
yobs = A[:, 1:end] # 观测值, 
# inds_obs = [2:7; 9; 11] # 观测SM对应的层
set_option!(; yobs, θ_surf, ibeg, same_layer=false)
options

soil = init_soil(; θ0, soil_type=7)
lower, upper = SM_paramBound(soil)
theta0 = SM_param2theta(soil)
ysim = model_sim(theta0;)

goal(theta0; ibeg)
plot_result(theta0)

begin
  # plot_result(theta0)
  f(theta) = goal(theta; ibeg)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=Int(2e4))
end

# serialize("theta", theta)
theta = deserialize("$(@__DIR__)/theta")
SM_UpdateParam!(soil, theta)

plot_result(theta)
## 最后，提出一种参数化方案
