using Plots, Printf
gr(framestyle=:box)


function plot_sim(i; ibeg, ysim=nothing)
  t = d[:, :time]
  depth = round(soil.z[ibeg+i-1] * 100, digits=0)

  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, format.(ticks, "mm-dd")

  p = plot(title="Layer Depth: $depth cm"; xticks,
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom)
  plot!(p, t, yobs[:, i], label="OBS")
  plot!(p, t, ysim[:, i], label="SIM")

  return p
end


function plot_result(theta)
  ysim = model_sim(theta)
  n = size(ysim, 2)
  plot([plot_sim(i; ibeg, ysim) for i in 1:n]..., size=(1200, 600))
end
