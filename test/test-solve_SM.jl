using SoilDifferentialEquations, OrdinaryDiffEq, Test


function data_loader_soil(; dt=60)
  N = 150
  _param = (θ_sat=0.287, θ_res=0.075, Ksat=34 / 3600, α=0.027, n=3.96, m=1.0)
  param = Init_SoilWaterParam(N, _param...; use_m=true)
  param_water = ParamVanGenuchten(; _param...)

  Δz = fill(0.01, N)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  θ = fill(0.1, N)
  ψ = Retention_ψ.(θ; par=param_water)
  θ0 = 0.267
  ψ0 = Retention_ψ(θ0; par=param_water)

  sink = ones(N) * 0.3 / 86400 # [cm s⁻¹], 3mm/d, 蒸发速率
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, θ, ψ, θ0, ψ0, dt, sink, param_water, param)
  return soil
end

begin
  # 4 hours
  dt = 60 * 6
  soil = data_loader_soil(; dt)
  ntime = round(Int, 3600 * 4 / dt)
  θ_surf = fill(0.267, ntime)
  ysim_bonan = solve_SM_Bonan(soil, θ_surf)

  soil = data_loader_soil(; dt)
  ysim_ode = solve_SM_ODE(soil, θ_surf; solver=Tsit5())

  @test maximum(abs.(ysim_bonan[end, :] - ysim_ode[end, :])) <= 0.03
end

begin
  # 4 hours
  dt = 60 * 6
  soil = data_loader_soil(; dt)
  ntime = round(Int, 3600 * 4 / dt)
  θ_surf = fill(0.267, ntime)
  ysim_bonan = ModSim_SM(soil, θ_surf; method="Bonan")

  soil = data_loader_soil(; dt)
  ysim_ode = ModSim_SM(soil, θ_surf; method="Bonan", solver=Tsit5())
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
