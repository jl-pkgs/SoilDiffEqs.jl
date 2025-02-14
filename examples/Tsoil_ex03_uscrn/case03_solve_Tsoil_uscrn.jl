using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using OrdinaryDiffEq, Ipaper
import HydroTools: sceua, GOF, of_KGE, of_NSE
using Artifacts
includet("main_Tsoil.jl")


begin
  # :SOLARAD
  vars_TS = [:SUR_TEMP, :SOIL_TEMP_5, :SOIL_TEMP_10, :SOIL_TEMP_20, :SOIL_TEMP_50, :SOIL_TEMP_100]
  vars = [:T_HR_AVG; vars_TS...]

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)
end


function Tsoil_calib()
  # kw_solver = (; method="ODE", solver)
  kw_solver = (; method="Bonan")

  info = OrderedDict()
  for i in eachindex(sites)
    printstyled("[i = $i] site = $(sites[i]) \n", color=:green, bold=true)

    SITE = sites[i]
    d = df[df.site.==SITE, :][1:24*14, [:time; vars]]
    # t = DateTime.(d.time, "yyyy-mm-ddTHH:MM:SSZ")

    yobs_full = d[:, vars_TS] |> Matrix |> drop_missing
    ibeg = 2
    yobs = yobs_full[:, ibeg:end]

    Tsoil0 = yobs_full[1, :]
    Tsurf = yobs_full[:, ibeg]

    soil = init_soil(; Tsoil0, soil_type=7)
    x0 = [soil.κ; soil.cv]
    nlayer = length(soil.κ)
    lower = [fill(0.1, nlayer); fill(0.1, nlayer) * 1e6]
    upper = [fill(10.0, nlayer); fill(5.0, nlayer) * 1e6]

    f(theta) = goal(theta; yobs, Tsoil0, Tsurf, ibeg, kw_solver...)
    @time theta, feval, exitflag = sceua(f, x0, lower, upper; maxn=Int(5e4))

    soil = init_soil(; Tsoil0, soil_type=7)
    ysim = model_Tsoil_sim(soil, Tsurf, theta; kw_solver...)
    info[SITE] = -feval
  end
  info2 = DataFrame(; site=collect(keys(info)), NSE=collect(values(info)))
  return info2
end

info2 = Tsoil_calib()
# 50/54个模型，NSE表现超过0.9
info2[info2.NSE.<0.9, :]


# 一轮之后，温度发生了变化。
# soil = init_soil(; soil_type=7, ibeg)
# theta = [soil.κ; soil.cv]
# ysim = model_Tsoil_sim(soil, Tsurf, theta; )
# _n = length(soil.inds_obs)
# plot([plot_soil(i; ibeg) for i in 1:_n]..., size=(1200, 800))

# r = optimize(f, lower, upper, x0, Fminbox(inner_optimizer), options)
# theta = r.minimizer
