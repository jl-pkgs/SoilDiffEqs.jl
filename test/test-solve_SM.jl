## TODO: Bonan求解不稳定，会来回跳动
using SoilDifferentialEquations, OrdinaryDiffEqTsit5, Test
# using Plots

function init_soil(; dt=60)
  N = 150
  par = VanGenuchten(θ_sat=0.287, θ_res=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1.0)
  param = SoilParam(N, par; use_m=true)

  Δz = fill(0.01, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.1, N)
  ψ = Retention_ψ.(θ; par)
  θ0 = 0.267
  ψ0 = Retention_ψ(θ0; par)

  sink = ones(N) * 0.3 / 86400 # [cm s⁻¹], 3mm/d, 蒸发速率
  sink = zeros(N) * 0.3 / 86400 # [cm s⁻¹], 3mm/d, 蒸发速率
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, θ, ψ, θ0, ψ0, dt, sink, param)
  return soil
end

begin
  # 4 hours
  dt = 60 * 6
  soil = init_soil(; dt)
  ntime = round(Int, 3600 * 4 / dt)
  θ_surf = fill(0.267, ntime)
  ysim_bonan = solve_SM_Bonan(soil, θ_surf)

  soil = init_soil(; dt)
  ysim_ode = solve_SM_ODE(soil, θ_surf; solver=Tsit5())

  @test maximum(abs.(ysim_bonan[end, :] - ysim_ode[end, :])) <= 0.03
end

# begin
#   gr(framestyle=:box)
#   t_max = 3600 * 4
#   t = dt:dt:t_max
#   x = 1800:1800:t_max
#   inds = indexin(x, t)
  
#   ps = []
#   for i in inds
    
#     p = plot(; title="Time = $(i * dt / 3600) h")
#     y_bonan = ysim_bonan[i, :]
#     y_ode = ysim_ode[i, :]
#     z = soil.z₊ₕ
    
#     plot!(p, y_bonan, z, label="Bonan")
#     plot!(p, y_ode, z, label="ODE")
#     push!(ps, p)
#   end
#   plot(ps..., size=(1200, 800), ylabel="Depth (m)", xlabel="θ")
# end

begin
  # 4 hours
  dt = 60 * 6
  soil = init_soil(; dt)
  ntime = round(Int, 3600 * 4 / dt)
  θ_surf = fill(0.267, ntime)
  ysim_bonan = ModSim_SM(soil, θ_surf; method="Bonan")

  soil = init_soil(; dt)
  ysim_ode = ModSim_SM(soil, θ_surf; method="ODE", solver=Tsit5())
  @test maximum(abs.(ysim_bonan[end, :] - ysim_ode[end, :])) <= 0.03
end

# begin
#   function plot_sm(i)
#     title = "Layer $i"
#     plot(; title)    
#     plot!(ysim_bonan[:, i]; label = "Bonan")
#     plot!(ysim_ode[:, i]; label = "ODE")
#   end
#   layers = [1, 5, 10, 20, 50, 100, 150]
#   plot([plot_sm(i) for i in layers]..., size=(1200, 800))
# end
