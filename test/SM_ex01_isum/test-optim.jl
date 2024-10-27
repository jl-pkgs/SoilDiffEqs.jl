# 优化模型参数
using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using OrdinaryDiffEq
import HydroTools: sceua, GOF, of_KGE, of_NSE
import NetCDFTools: approx
using Optim

includet("main_Tsoil.jl")

# Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
begin
  soil = init_soil(; soil_type=7)
  soil.Tsoil = 
  # x0 = [soil.κ[1]; soil.cv[1]]
  # lower = [0.1; 0.1 * 1e6]
  # upper = [10.0; 5.0 * 1e6]
  x0 = [soil.κ; soil.cv]
  _n = soil.n
  lower = [fill(0.1, _n); fill(0.1, _n) * 1e6]
  upper = [fill(10.0, _n); fill(5.0, _n) * 1e6]

  f(theta) = goal(theta; ibeg)
  # inner_optimizer = GradientDescent()
  # options = Optim.Options(show_trace=true)
  # r = optimize(f, lower, upper, x0, Fminbox(inner_optimizer), options)
  @time theta, feval, exitflag = sceua(f, x0, lower, upper; maxn=Int(1e5))
end

soil = init_soil(; soil_type=7)
x0 = [soil.κ; soil.cv]
# theta = r.minimizer
ysim = model_sim(x0; ibeg)

plot(
  [plot_soil(i) for i in 1:8]...,
  size=(1200, 800),
)
# of_NSE(yobs, ysim)

# z, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
# nlayer = length(dz)
# soil = Soil(dz)
