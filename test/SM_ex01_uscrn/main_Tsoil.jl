using Plots, Printf
gr(framestyle=:box)


function plot_sim(i; ibeg=2)
  i2 = i + ibeg - 1
  title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)
  plot(; title)
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

# inner_optimizer = GradientDescent()
# options = Optim.Options(show_trace=true)
