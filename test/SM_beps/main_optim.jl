using SoilDifferentialEquations
using SoilDifferentialEquations.GlobalOptions

# solver = Tsit5()
# solver = Rosenbrock23()
# solver = Rodas5(autodiff=false)
z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100# 第一层是虚拟的

## TODO: 要想一种重复利用数据的方法
# used a global variable: `options`
# soil = init_soil()
function init_soil(; θ0, dt=3600.0, soil_type=7)
  (; method_retention, same_layer, ibeg) = options
  # dz = [2.5, 5, 5, 15, 45, 55]
  z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
  Δz = cal_Δz(z)
  N = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.2, N)
  θ[ibeg:end] .= θ0
  param_water = get_soilpar(soil_type)
  param = Init_SoilWaterParam(N, Vector(param_water)...;
    use_m=false, method=method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, param, param_water)
  soil.param.θ_sat .= 0.30
  soil.param.θ_res .= 0.03
  soil
end


function model_sim(theta)
  (; θ_surf, θ0) = options
  soil = init_soil(; θ0, soil_type=8)
  SM_UpdateParam!(soil, theta) # update param

  ysim = soil_moisture_BEPS(soil, θ0, θ_surf)
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
      ∑ += - GOF(obs, sim).NSE
  end, ibeg:ncol)
  ∑ / n # mean of NSE
end
