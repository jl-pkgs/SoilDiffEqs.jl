using Plots, Dates

function plot_sm(band)  
  t = d[:, :time]
  x = d[:, band]

  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, format.(ticks, "mm-dd")

  p = plot(title=string(band); xticks, 
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom,
  )
  plot!(p, t, x)
  return p
end

Figure1_raw = plot([plot_sm(band) for band in vars_SM]..., size=(1000, 600))
# findall(isnan.(yobs_full))

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
