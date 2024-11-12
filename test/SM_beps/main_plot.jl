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


function plot_result(ysim)
  # ysim = model_sim(theta)
  plot([plot_sim(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
end
