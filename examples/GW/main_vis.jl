function plot_soil_profile(soil)
  z = soil.z[1:N]
  p = plot(; size=(1000, 600), legendposition=:bottomleft,
    xlabel="θ (m³ m⁻³)", ylabel="Depth (m)", title="(a) Soil Profile (zwt = $zwt m)")

  days = [1, 2, 3, 5, 7, 20, 30, 60, 90, 120, 150, 180, 210]
  for day in days
    t = day * 24
    plot!(p, SM[t, :], z, label="day=$day")
  end
  plot!(p, soil.θE[1:N], z, label="θE", color=:black, linewidth=2)
  p
end


## 不同深度的时间序列
function plot_depth_timeseries(soil)
  depths = [2.5, 10, 20, 30, 40, 50, 70, 100, 130, 180]
  _z = round.(-soil.z[1:N] * 100, digits=1)
  layers = indexin(depths, _z)
  # layers = round.(Int, depths / (Δ*100))

  p = plot(; title="(b) θ time-series (zwt = $zwt m)", xlabel="Date", ylabel="θ (m³ m⁻³)")
  t = 1:24*90
  t = 1:length(dates)
  for i in eachindex(depths)
    l = layers[i]
    depth = depths[i]
    # depth = round(Int, l * Δ * 100)
    plot!(p, dates[t], SM[t, l], label="$(depth)cm")
  end
  p
end


function plot_θ(soil)
  N = soil.N
  z = soil.z[1:N]
  zwt = soil.zwt

  plot(title="zwt = $zwt m", xlabel="θ (m³ m⁻³)", ylabel="z [m]", legend=:topright)
  plot!(θ, z, label="θ_init")
  plot!(soil.θ, z, label="θ_next")
end
