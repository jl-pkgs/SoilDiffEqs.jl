using SoilDifferentialEquations, Test
import RTableTools: fread

const Z_OBS = [0.05, 0.10, 0.20, 0.50, 1.00]
const DT = 3600.0

# data
raw = fread(f_SM_Batesville)
yobs_raw = Matrix(raw[:, 3:end])

layers = get_clm5_layers()
z_clm = abs.(layers.z[1:length(DZ_CLM5)])

yobs_interp = interp_obs_to_layer(yobs_raw, Z_OBS, z_clm)
θ_surf = yobs_raw[:, 1]
θ0_initial = yobs_interp[1, :]

set_option!(; yobs=yobs_interp, θ_surf, ibeg=1, same_layer=true)

function init_soil_clm5(; θ0=nothing, dt=DT, soil_type=7, use_m=false)
  (; method_retention, same_layer) = options
  layers = get_clm5_layers()
  (; z, z₋ₕ, z₊ₕ, Δz, Δz₊ₕ) = layers
  N = length(Δz)

  θ = fill(0.2, N)
  θ0 !== nothing && (θ .= θ0)

  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par; use_m, method_retention, same_layer)

  Soil{Float64}(; N, ibeg=1, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
end

function model_sim(theta)
  soil = init_soil_clm5(; θ0=θ0_initial, soil_type=8)
  SM_UpdateParam!(soil, theta)

  method = options.method_solve
  if method == "Bonan"
    solve_SM_Bonan(soil, θ_surf)
  else
    solve_SM_ODE(soil, θ_surf; solver=Tsit5())
  end
end

function goal(theta)
  yobs = options.yobs
  ysim = model_sim(theta)

  ncol = size(yobs, 2)
  loss = 0.0
  for i in 1:ncol
    obs = @view yobs[:, i]
    sim = @view ysim[:, i]
    loss += -GOF(obs, sim).NSE
  end
  loss / ncol
end

println("Running simulation with CLM5 layers...")
θ_test = SM_param2theta(init_soil_clm5(; θ0=θ0_initial))
score = goal(θ_test)
println("NSE Score (negative mean): $score")

@testset "CLM5 Simulation" begin
  @test !any(isnan.(θ_test))
  @test score isa Real
  @test score < 100.0
end
