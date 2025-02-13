using SoilDifferentialEquations, Test, Dates
using SoilDifferentialEquations.GlobalOptions
import RTableTools: fread
include("main_optim.jl")


begin
  d = fread(f_SM_Batesville)
  ibeg = 2
  yobs_full = d[:, 3:end] |> Matrix #|> drop_missing
  yobs = yobs_full[:, max(ibeg - 1, 1):end]
  θ0 = yobs_full[1, max(ibeg - 1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  set_option!(; yobs, θ_surf, ibeg, same_layer=true)
  options
end


function test_ModSim(; method_retention, maxn=2_000, kw...)
  set_option!(; method_retention, kw...)

  # [5, 10, 20, 50, 100]
  soil = init_soil(; θ0, soil_type=7)
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  goal(theta0)
  # ysim = model_sim(theta0)

  f(theta) = goal(theta; ibeg=1)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn)
  -feval
end


# @testset "ModSim_SM" 
@profview begin
  printstyled("[ModSim_SM]: Campbell ... \n", color=:blue, bold=true)
  test_ModSim(; method_retention="Campbell", same_layer=false) >= 0.32
  # printstyled("[ModSim_SM]: van Genuchten ... \n", color=:blue, bold=true)
  
  # @test test_ModSim(; method_retention="van_Genuchten", same_layer=false) >= 0.10
  # @test test_ModSim(; method_retention="Campbell", same_layer=false, maxn=2000) >= 0.60
  # @test test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=5000) >= 0.65
end
# @profview test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=5000)
