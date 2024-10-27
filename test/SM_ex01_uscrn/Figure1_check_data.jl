Figure1_raw = plot([plot_SM(band) for band in vars_SM]..., size=(1000, 600))

begin
  soil = init_soil(; soil_type=7)
  method = "ODE"
  method = "Bonan"

  if method == "Bonan"
    ysim = solve_Tsoil_Bonan(soil, TS0)
  elseif method == "ODE"
    solver = Tsit5()
    # solver = Rosenbrock23()
    # solver = Rodas5(autodiff=false)
    ysim = solve_Tsoil_ODE(soil, TS0; solver)
  end
  # savefig("case01_Tsoil_CUG_$(method).png")
end
