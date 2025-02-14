using SoilDifferentialEquations, Test, Dates
import RTableTools: fread
include("main_optim.jl")

begin
  d = fread(f_SM_Batesville)
  ibeg = 2
  yobs_full = d[:, 3:end] |> Matrix #|> drop_missing
  yobs = yobs_full[:, max(ibeg - 1, 1):end] # [time, depth]
  θ0 = yobs_full[1, max(ibeg - 1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  set_option!(; yobs, θ_surf, ibeg, same_layer=true)
  options
end


function test_ModSim(; method_retention, maxn=4_000, use_m=false, kw...)
  set_option!(; method_retention, kw...)
  # [5, 10, 20, 50, 100]
  soil = init_soil(; θ0, soil_type=7, use_m)
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  npar = length(theta0)
  println("npar = $npar")

  goal(theta0)
  # ysim = model_sim(theta0)
  f(theta) = goal(theta; ibeg=1)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn)
  -feval
end

@testset "ModSim_SM" begin
  @test test_ModSim(; method_retention="Campbell", same_layer=false) >= 0.7
  @test test_ModSim(; method_retention="van_Genuchten", same_layer=false) >= 0.1
end
