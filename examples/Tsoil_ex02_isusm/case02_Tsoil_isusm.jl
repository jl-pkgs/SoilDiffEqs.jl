using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using OrdinaryDiffEq
import HydroTools: sceua, GOF, of_KGE, of_NSE
import NetCDFTools: approx
import Ipaper: set_seed

includet("main_Tsoil.jl")

solver = Tsit5()
# solver = Rosenbrock23()
# solver = Rodas5(autodiff=false)

## 加载数据
begin
  d = fread("$dir_soil/data/isusm_TS_202207.csv")[1:24*10, :]
  t = d.time

  ibeg = 3
  yobs_full = Matrix(d[:, 4:end]) # 共有12层土壤数据

  set_seed(1)
  Tsurf = yobs_full[:, ibeg] # 从第10层开始
  Tsoil0 = approx(inds_obs, yobs_full[1, :], 1:nlayer)

  soil = init_soil(; soil_type=7)
  theta0 = [soil.κ; soil.cv]
  # theta = theta0

  yobs = yobs_full[:, ibeg:end]
  ysim = model_Tsoil_sim(soil, Tsurf, theta; method="ODE", solver)
  # ysim = model_Tsoil_sim(soil, Tsurf, theta0; method="Bonan")
end

## 结论是初始值很重要
begin
  plot(
    [plot_soil(i) for i in 1:10]...,
    size=(1200, 800))
  # savefig("case02_Bonan_raw.png")
end


## 优化模型参数
begin
  x0 = [soil.κ; soil.cv]
  lower = [fill(0.1, nlayer); fill(1.0, nlayer) * 1e6]
  upper = [fill(10.0, nlayer); fill(5.0, nlayer) * 1e6]

  f(theta) = goal(theta; method="ODE", solver)
  # f(theta) = goal(theta; method="Bonan")
  @time theta, feval, exitflag = sceua(f, x0, lower, upper; maxn=Int(5e4))
end
