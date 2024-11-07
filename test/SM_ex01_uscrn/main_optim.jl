using Plots, Printf
gr(framestyle=:box)

function plot_sim(i; ibeg, ysim=nothing)
  band = vars_SM[i]
  t = d[:, :time]
  x = d[:, band]
  k = i - ibeg + 1

  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, format.(ticks, "mm-dd")

  p = plot(title=string(band); xticks,
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom)
  plot!(p, t, x, label="OBS")
  if k >= 1 && ysim !== nothing
    # i2 = i - 1
    # title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)
    plot!(p, t, ysim[:, k], label="SIM")
  end
  return p
end


z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100# 第一层是虚拟的

function init_soil(; θ0, dt=3600.0, ibeg=2, soil_type=7, same_layer=true)
  # dz = [2.5, 5, 5, 15, 45, 55]
  z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
  Δz = cal_Δz(z)
  N = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  # m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  θ = fill(0.2, N)
  θ[ibeg:end] .= θ0
  param_water = get_soilpar(soil_type)
  param = Init_SoilWaterParam(N, Vector(param_water)...;
    use_m=false,
    method="Campbell",
    # method="van_Genuchten",
    same_layer)
  Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, param, param_water)
end


function model_sim(theta; method="Bonan", same_layer=true)
  soil = init_soil(; θ0, soil_type=8, ibeg, same_layer)
  SM_UpdateParam!(soil, theta) # update param

  if method == "Bonan"
    ysim = solve_SM_Bonan(soil, θ_surf)
  else
    # solver = Rosenbrock23()
    # solver = Rodas5(autodiff=false)
    solver = Tsit5()
    ysim = solve_SM_ODE(soil, θ_surf; solver)
  end
  return ysim
end

function goal(theta; method="Bonan", same_layer=true, ibeg=1)
  ysim = model_sim(theta; method, same_layer)
  obs = yobs[:, ibeg:end]
  sim = ysim[:, ibeg:end]
  -GOF(obs[:], sim[:]).NSE
end
