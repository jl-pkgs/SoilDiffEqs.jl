using Plots

function plot_tem(band)
  x = d[:, band] * 1.0
  p = plot(title=string(band))
  plot!(p, t, x)
  return p
end

Figure1_raw = plot([plot_tem(band) for band in vars]..., size=(1000, 600))
Figure1_raw

# findall(isnan.(yobs_full))

begin
  soil = init_soil(; soil_type=7)
  method = "ODE"
  method = "Bonan"

  if method == "Bonan"
    ysim = solve_Tsoil_Bonan(soil, Tsurf)
  elseif method == "ODE"
    solver = Tsit5()
    # solver = Rosenbrock23()
    # solver = Rodas5(autodiff=false)
    ysim = solve_Tsoil_ODE(soil, Tsurf; solver)
  end
  # savefig("case01_Tsoil_CUG_$(method).png")
end
