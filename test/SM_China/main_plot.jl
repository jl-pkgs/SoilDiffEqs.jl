using Plots, Printf, Dates
gr(framestyle=:box)


function plot_sim(i; ysim=nothing, ignore...)
  t = dates
  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Month(1):time_max
  xticks = ticks, format.(ticks, "mm-dd")

  title = @sprintf("%dcm", depths[i])
  p = plot(; title, 
    xticks,
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom)

  plot!(p, t, yobs[:, i], label="OBS")
  plot!(p, t, ysim[:, i], label="SIM")
  return p
end


function plot_result(theta)
  ysim = model_sim(theta)
  nlys = size(ysim, 2)
  plot([plot_sim(i; ysim) for i in 1:nlys]..., size=(1400, 800))
end

# plot_result(theta)
## 部分土壤层，观测存在错误

# begin
#   plot_sim(1; ysim)
#   plot_result(theta)
# end
