using SoilDifferentialEquations, Test, Dates
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

# test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=4_000, use_m=false)
# # npar = 30
# # Iteration =   0, nEvals = 305, Best Cost = 1.08955
# # Iteration =   1, nEvals = 724, Best Cost = 1.08955
# # Iteration =   2, nEvals = 1205, Best Cost = 0.05273
# # Iteration =   3, nEvals = 1727, Best Cost = 0.05273
# # Iteration =   4, nEvals = 2253, Best Cost = -0.10765
# # Iteration =   5, nEvals = 2800, Best Cost = -0.33195
# # Iteration =   6, nEvals = 3351, Best Cost = -0.33195
# # Iteration =   7, nEvals = 3891, Best Cost = -0.33195
# # Iteration =   8, nEvals = 4434, Best Cost = -0.37321
# #   21.189072 seconds (3.77 M allocations: 5.693 GiB, 3.15% gc time)

# test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=4_000, use_m=true)
# # Iteration =   0, nEvals = 365, Best Cost = -0.01120
# # Iteration =   1, nEvals = 863, Best Cost = -0.01120
# # Iteration =   2, nEvals = 1457, Best Cost = -0.01120
# # Iteration =   3, nEvals = 2054, Best Cost = -0.01120
# # Iteration =   4, nEvals = 2666, Best Cost = -0.11453
# # Iteration =   5, nEvals = 3303, Best Cost = -0.19156
# # Iteration =   6, nEvals = 3937, Best Cost = -0.31047
# # Iteration =   7, nEvals = 4570, Best Cost = -0.47214

@testset "ModSim_SM" begin
  @test test_ModSim(; method_retention="Campbell", same_layer=false) >= 0.7
  @test test_ModSim(; method_retention="van_Genuchten", same_layer=false) >= 0.1
end
