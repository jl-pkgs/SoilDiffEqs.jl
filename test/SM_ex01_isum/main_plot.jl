using Plots
gr(framestyle=:box)


function plot_soil(i)
  plot(title="layer $i")
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

function plot_obs(i)
  plot(title="layer $i")
  plot!(t, yobs[:, i], label="OBS")
  # plot!(t, TS_sim[:, i], label="SIM")
end
