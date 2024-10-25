# 优化模型参数
using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using DifferentialEquations
import HydroTools: sceua, GOF, of_KGE, of_NSE
using Optim
includet("main_Tsoil.jl")

d = fread("$dir_soil/data/CUG_TS_202306.csv")

t = d.time
A = Matrix(d[:, 2:end]) #|> drop_missing 

ibeg = 2
TS0 = A[:, ibeg]
yobs = A[:, ibeg:end]

begin
  soil = init_soil(; soil_type=7)
  # x0 = [soil.κ[1]; soil.cv[1]]
  # lower = [0.1; 0.1 * 1e6]
  # upper = [10.0; 5.0 * 1e6]
  x0 = [soil.κ; soil.cv]
  lower = [fill(0.1, 9); fill(0.1, 9) * 1e6]
  upper = [fill(10.0, 9); fill(5.0, 9) * 1e6]

  f(theta) = goal(theta; ibeg=2)
  
  # inner_optimizer = GradientDescent()
  # options = Optim.Options(show_trace=true)
  # r = optimize(f, lower, upper, x0, Fminbox(inner_optimizer), options)
  @time theta, feval, exitflag = sceua(f, x0, lower, upper; maxn=Int(1e5))
end

# theta = r.minimizer
ysim = model_sim(theta; ibeg)

plot(
  [plot_soil(i) for i in 1:8]...,
  size=(1200, 800),
)
# of_NSE(yobs, ysim)

# z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
# nlayer = length(dz)
# soil = Soil(dz)
