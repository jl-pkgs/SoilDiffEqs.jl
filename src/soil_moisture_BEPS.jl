# ψ: 负
# z: 向下为负
function soil_moisture_BEPS(soil::Soil{FT}; method="ψ0") where {FT}
  kstep = soil.dt # 1 hour
  ∑t = 0.0

  # TODO: 输入边界条件
  θ = soil.θ

  while ∑t < kstep
    Q0 = soil_WaterFlux!(soil, θ; method) # update Q
    Qmax = maximum(abs.(soil.Q))
    Qmax = max(Qmax, abs(Q0))

    Δt = guess_step(Qmax) # this_step
    ∑t += Δt
    ∑t > kstep && (Δt -= (∑t - kstep))

    soil_Updateθ!(soil, Δt) # update θ during Δt
  end
end

"""
    guess_step(max_Fb)

    如果流速过快，则减小时间步长

# Arguments
- `max_Fb`: [cm s-1]

# Notes
[m s-1] -> 1000*[mm s-1] -> 1000*[kg m-2 s-1]
"""
function guess_step(max_Fb)
  # this constraint is too large
  # 需要调研，暴雨期间土壤的下渗速率，或者说`Ksat`
  if max_Fb > 1.0e-5 # 864 mm/day
    dt = 1.0
  elseif max_Fb > 1.0e-6 # 86.4 mm/day
    dt = 30.0 # seconds
  else
    dt = 360.0
  end
  dt
end

export soil_moisture_BEPS
