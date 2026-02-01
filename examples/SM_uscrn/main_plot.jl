using Plots, Printf, Dates
gr(framestyle=:box)

function plot_sim(i; ysim, yobs, dates, var_name)
  t = dates
  x = yobs[:, i]
  
  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, format.(ticks, "mm-dd")
  
  p = plot(title=var_name; xticks,
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom)
  plot!(p, t, x, label="OBS")
  if ysim !== nothing && i <= size(ysim, 2)
    plot!(p, t, ysim[:, i], label="SIM")
  end
  return p
end

function plot_result(; ysim, yobs, dates, var_names, filename=nothing)
  n = min(length(var_names), size(yobs, 2))
  plots = []
  for i in 1:n
    # 只绘制有模拟结果的层
    if i <= size(ysim, 2)
      push!(plots, plot_sim(i; ysim, yobs, dates, var_name=var_names[i]))
    else
      push!(plots, plot_sim(i; ysim=nothing, yobs, dates, var_name=var_names[i]))
    end
  end
  p = plot(plots..., size=(1200, 600))
  
  if filename !== nothing
    savefig(p, filename)
    println("Plot saved: $filename")
  end
  p
end
