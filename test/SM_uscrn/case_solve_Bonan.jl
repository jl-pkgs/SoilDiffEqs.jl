using SoilDifferentialEquations, Test, Dates
import RTableTools: fread

begin
  d = fread(f_SM_Batesville)
  ibeg = 3
  yobs_full = d[:, 3:end] |> Matrix #|> drop_missing

  yobs = yobs_full[:, max(ibeg - 1, 1):end] # [time, depth]
  θ0 = yobs_full[1, max(ibeg - 1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  set_option!(; yobs, θ_surf, ibeg, same_layer=true)
  options
end


function init_soil(; θ0=0.3, dt=3600.0, soil_type=7, use_m=false)
  (; method_retention, same_layer, ibeg) = options
  # dz = [2.5, 5, 5, 15, 45, 55]
  z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
  Δz = cal_Δz(z)
  N = length(Δz)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.2, N)
  θ[ibeg:end] .= θ0
  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par;
    use_m, method_retention, same_layer)
  Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
end
# soil = init_soil()

function model_sim(theta)
  soil = init_soil(; θ0, soil_type=8)
  SM_UpdateParam!(soil, theta) # update param

  θ_surf = options.θ_surf
  method = options.method_solve
  if method == "Bonan"
    ysim = solve_SM_Bonan(soil, θ_surf)
  else
    solver = Tsit5()
    ysim = solve_SM_ODE(soil, θ_surf; solver)
  end
  return ysim
end

function goal(theta; ibeg=1)
  yobs = options.yobs
  ysim = model_sim(theta)

  ncol = size(yobs, 2)
  n = ncol - ibeg + 1
  ∑ = 0.0
  map(i -> begin
      obs = @view yobs[:, i]
      sim = @view ysim[:, i]
      ∑ += -GOF(obs, sim).NSE
    end, ibeg:ncol)
  ∑ / n # mean of NSE
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

test_ModSim(; method_retention="van_Genuchten", same_layer=false, maxn=20_000)


test_ModSim(; method_retention="Campbell", same_layer=false, maxn=5_000)

