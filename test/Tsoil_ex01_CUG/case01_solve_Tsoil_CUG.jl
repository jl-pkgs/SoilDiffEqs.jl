# TODO: 优化模型参数
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
Tsoil0 = A[1, :]

soil = init_soil(; soil_type=7)
# soil.κ[5:9] .= 4.0
ysim = solve_Tsoil_ODE(soil, TS0; ibeg)

plot(
  [plot_soil(i; ibeg) for i in 1:8]...,
  size=(1200, 800),
)

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
  @time theta, feval, exitflag = sceua(f, x0, lower, upper; maxn=Int(1 * 1e3))
end
ysim = model_sim(theta; ibeg)

# theta0 = [fill(0.2, 9); fill(0.1, 9) * 1e6]
# goal(theta; ibeg=2)
# theta = r.minimizer
