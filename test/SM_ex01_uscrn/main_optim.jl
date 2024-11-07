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

# θ_sat, θ_res, Ksat, α, n, m
LOWER = [0.25, 0.03, 0.002 / 3600, 0.002, 1.05, 0.1]
UPPER = [0.50, 0.20, 60.0 / 3600, 0.300, 4.00, 10.0]

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
  UpdateSoilParam!(soil, theta) # update param

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


function UpdateSoilParam!(soil::Soil{T}, theta::AbstractVector{T}) where {T<:Real}
  if soil.param.same_layer
    soil.param.θ_sat .= theta[1]
    soil.param.θ_res .= theta[2]
    soil.param.Ksat .= theta[3]
    soil.param.α .= theta[4]
    soil.param.n .= theta[5]
    soil.param.m .= theta[6]
  else
    N = soil.N
    soil.param.θ_sat .= theta[1:N]
    soil.param.θ_res .= theta[N+1:2N]
    soil.param.Ksat .= theta[2N+1:3N]
    soil.param.α .= theta[3N+1:4N]
    soil.param.n .= theta[4N+1:5N]
    soil.param.m .= theta[5N+1:6N]
  end
  return nothing
end

function param2theta(soil)
  (; θ_sat, θ_res, Ksat, α, n, m) = soil.param
  if soil.param.same_layer
    return [θ_sat[1], θ_res[1], Ksat[1], α[1], n[1], m[1]]
  else
    return [θ_sat; θ_res; Ksat; α; n; m]
  end
end

function get_bound(soil)
  if soil.param.same_layer
    return LOWER, UPPER
  else
    return repeat(LOWER; inner=soil.N), repeat(UPPER; inner=soil.N)
  end
end
