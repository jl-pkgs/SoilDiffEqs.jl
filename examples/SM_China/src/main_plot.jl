using Plots, Dates

# ========== Plot Functions ==========
# 所有全局变量转为函数参数
# yobs 和 ysim 从第一层开始对应，depths 从 ibeg 开始对应
plot_sim(i; ysim, yobs, dates, depths) = begin
  depth = round(depths[i]; digits=0)
  t = dates
  time_min, time_max = extrema(t)
  ticks = time_min:Month(1):time_max
  
  p = plot(t, yobs[:, i]; label="OBS", title="$(Int(depth))cm",
    xticks=(ticks, format.(ticks, "mm-dd")), xrot=30,
    framestyle=:box, tickfonthalign=:center)
  plot!(p, t, ysim[:, i]; label="SIM")
end

plot_result(; ysim, yobs, dates, depths, ibeg, filename=nothing) = begin
  nlys = size(ysim, 2)
  # yobs 和 ysim 从第1层开始，对应 depths[ibeg-1:end]
  plot_depths = depths[ibeg-1:ibeg-1+nlys-1]
  plots = [plot_sim(i; ysim, yobs, dates, depths=plot_depths) for i in 1:nlys]
  p = plot(plots..., size=(1400, 800))
  
  if filename !== nothing
    savefig(p, filename)
    println("Plot saved: $filename")
  end
  p
end
