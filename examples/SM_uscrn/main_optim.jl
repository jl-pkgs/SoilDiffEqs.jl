using SoilDifferentialEquations
z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100# 第一层是虚拟的

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

# soil_type = 7
# method_retention = "van_Genuchten"
# same_layer = false
# N = 10
# par = get_soilpar(soil_type; method_retention)
# param = SoilParam(N, par;
#   use_m=false, method_retention, same_layer)

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
      obs = yobs[:, i]
      sim = ysim[:, i]
      ∑ += -GOF(obs, sim).NSE
    end, ibeg:ncol)
  ∑ / n # mean of NSE
end
