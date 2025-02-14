using SoilDifferentialEquations
using SoilDifferentialEquations.GlobalOptions

# solver = Tsit5()
# solver = Rosenbrock23()
# solver = Rodas5(autodiff=false)
# z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100# 第一层是虚拟的

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


function goal(theta;)
  yobs = options.yobs
  ysim = model_sim(theta)

  ncol = size(yobs, 2)
  n = ncol - ibeg + 1
  ∑ = 0.0
  
  for i in 1:ncol
    obs = @view yobs[:, i]
    sim = @view ysim[:, i]
    ∑ += -of_KGE(obs, sim)
  end
  return ∑ / n # mean of NSE
end
