using SoilDifferentialEquations, Test, Dates
using SoilDifferentialEquations.GlobalOptions
import RTableTools: fread
include("main_optim.jl")

GlobalOptions.options = Options()
options = GlobalOptions.options

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


function test_ModSim(; method_retention, maxn=5_000, kw...)
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


@testset "ModSim_SM" begin
  @test test_ModSim(; method_retention="Campbell", same_layer=false) >= 0.32
  @test test_ModSim(; method_retention="van_Genuchten", same_layer=false) >= 0.10
  # @test test_ModSim(; method_retention="Campbell", same_layer=false, maxn=2000) >= 0.60
  # @test test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=5000) >= 0.65
end
# @profview test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=5000)


# Iteration =   0, nEvals = 245, Best Cost = -0.23514
# Iteration =   1, nEvals = 598, Best Cost = -0.49180
# Iteration =   2, nEvals = 1004, Best Cost = -0.52860
# Iteration =   3, nEvals = 1411, Best Cost = -0.55392
# Iteration =   4, nEvals = 1804, Best Cost = -0.59312
# Iteration =   5, nEvals = 2167, Best Cost = -0.67321

# Iteration =   0, nEvals = 365, Best Cost = 2.70014
# Iteration =   1, nEvals = 869, Best Cost = 1.41141
# Iteration =   2, nEvals = 1479, Best Cost = 0.57537
# Iteration =   3, nEvals = 2134, Best Cost = -0.07305
# Iteration =   4, nEvals = 2746, Best Cost = -0.07305
# Iteration =   5, nEvals = 3339, Best Cost = -0.33710
# Iteration =   6, nEvals = 3894, Best Cost = -0.45565
# Iteration =   7, nEvals = 4439, Best Cost = -0.59242
# Iteration =   8, nEvals = 4982, Best Cost = -0.65706
# Iteration =   9, nEvals = 5561, Best Cost = -0.70967
