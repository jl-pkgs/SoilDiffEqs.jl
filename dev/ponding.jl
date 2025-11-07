# 蓄满产流
function ponding_simple(h_pond_current::T, P::T, inf::T, dt::T;
  h_pond_max::T=T(0.5), ignored...) where {T<:Real}

  h_pond = h_pond_current + (P - inf) * dt  # [cm]
  # 产流计算
  if h_pond > h_pond_max
    runoff = h_pond - h_pond_max
    h_pond_new = h_pond_max
  else
    runoff = 0.0
    h_pond_new = max(h_pond, 0.0)  # 确保非负
  end
  return runoff, h_pond_new
end


# 大部分变量都采用cm
function ponding_mixed(h_pond_current::T, P::T, inf::T, dt::T;
  h_pond_min::T=T(0.1), # 地表自由水的蓄水容量, [cm]
  h_pond_max::T=T(2.0), 
  S0::T=T(0.01),        # 坡面比降, [-]
  n_manning::T=T(0.15), # 曼宁粗糙度
  L::T=T(100.0),   # 网格长度, cm
  ignored...) where {T<:Real}

  h_pond = h_pond_current + (P - inf) * dt  # [cm], # 地表水量平衡

  # 分段处理
  if h_pond <= h_pond_min
    # 阶段1：微地形蓄水阶段，无径流
    runoff = 0.0
    h_pond_new = max(h_pond, 0.0)

  elseif h_pond <= h_pond_max
    # 阶段2：坡面汇流阶段，运动波方程
    h_excess = h_pond - h_pond_min
    h_m = h_excess / 100.0  # [cm] → [m]

    # 曼宁公式中众多参数，可以简化为一个参数
    v = (1.0 / n_manning) * (h_m^(2 / 3)) * sqrt(max(S0, 1e-6))  # [m/s]
    q = v * h_m  # [m²/s]
    q_cm_h = q * 3600.0 * 10000.0  # [cm²/h]
    drain_rate = q_cm_h / (L * 100)  # [cm/h], L [m] => [cm]

    # 排水量不能超过可用水量
    drainage = min(drain_rate * dt, h_excess)
    runoff = drainage
    h_pond_new = h_pond - drainage

  else
    # 阶段3：超蓄阶段，立即产流
    runoff = h_pond - h_pond_max
    h_pond_new = h_pond_max
  end

  return runoff, h_pond_new
end


@with_kw mutable struct ParamInfiltrationPhysics{FT}
  h_pond_min::T = 0.1 # 地表自由水的蓄水容量, [cm]，低于该数值不产流
  h_pond_max::T = 2.0 # 高于高数值全部产流
  n::FT = 0.15
  L::FT = 100.0 # [m]
  S0::FT = 0.01 # [m m-1]
end


@with_kw mutable struct ParamInfiltration{FT}
  h_pond_min::FT = 0.1 # 地表自由水的蓄水容量, [cm]，低于该数值不产流
  h_pond_max::FT = 2.0 # 高于高数值全部产流
  c_inf::FT = T
end

n = 0.15
L = 100.0
S0 = 0.01
c_inf = 1.0 / n * sqrt(S0) * 3600 * 10000.0 / (100L)
drainage = c_inf * h_excess ^ (5/3) # [1 cm h-1]
