using Plots, Printf, Dates
gr(framestyle=:box)

function plot_sim(i; ibeg, ysim, yobs, dates, depths)
  depth = round(depths[i]; digits=0)
  t = dates
  
  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, format.(ticks, "mm-dd")
  
  p = plot(title="Layer Depth: $depth cm"; xticks,
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom)
  plot!(p, t, yobs[:, i], label="OBS")
  plot!(p, t, ysim[:, i], label="SIM")
  
  return p
end

function plot_result(; ysim, yobs, dates, depths, filename=nothing)
  n = size(ysim, 2)
  p = plot([plot_sim(i; ibeg=1, ysim, yobs, dates, depths) for i in 1:n]..., size=(1200, 600))
  
  if filename !== nothing
    savefig(p, filename)
    println("Plot saved: $filename")
  end
  p
end
