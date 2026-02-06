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

plot_result(; ysim, yobs, dates, depths, fout=nothing) = begin
  nlys = size(ysim, 2)
  # yobs 和 ysim 从第1层开始，对应 depths[ibeg-1:end]
  plots = [plot_sim(i; ysim, yobs, dates, depths=depths) for i in 1:nlys]
  p = plot(plots..., size=(1400, 800))

  if fout !== nothing
    indir = dirname(fout)
    mkpath(indir)

    savefig(p, fout)
    println("Plot saved: $fout")
  end
  p
end

# 绘制原始观测数据（诊断用）
plot_obs(data_origin, dates, zs_obs; fout=nothing) = begin
  nlys = size(data_origin, 2)
  time_min, time_max = extrema(dates)
  ticks = time_min:Month(1):time_max

  plots = map(1:nlys) do i
    depth = round(zs_obs[i]; digits=0)
    plot(dates, data_origin[:, i]; label="OBS", title="$(Int(depth))cm",
      xticks=(ticks, format.(ticks, "mm-dd")), xrot=30, framestyle=:box)
  end
  p = plot(plots..., size=(1400, 800))

  !isnothing(fout) && (mkpath(dirname(fout)); savefig(p, fout); println("Plot saved: $fout"))
  p
end

