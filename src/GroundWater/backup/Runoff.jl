# https://github.com/CoLM-SYSU/CoLM202X/blob/master/main/MOD_Runoff.F90
module Runoff

export surface_runoff_simtop, subsurface_runoff_simtop, runoff_xinanjiang, runoff_simplevic

# =============================================================================
# 函数：surface_runoff_simtop
#
# 功能：根据 TOPMODEL 理论计算地表径流，返回一个三元组(rsur, rsur_se, rsur_ie)：
# - rsur    ：总地表径流（单位：mm H₂O/s）
# - rsur_se ：饱和超渗径流部分
# - rsur_ie ：入渗超渗径流部分
#
# 参数说明：
#   N    :: 土壤层数（整型）
#   wimp       :: 不渗透阈值（T 类型）
#   Ksat       :: 各层饱和时的水力传导率（mm H₂O/s），向量
#   Fmax       :: 饱和区域最大分数（无量纲，T 类型）
#   f          :: 饱和区域衰减因子（1/m，T 类型）
#   eff_porosity::各层有效孔隙率（T 类型），向量
#   icefrac    :: 各层冰的比例（无量纲），向量
#   zwt        :: 从土壤表面到水位的深度（m，T 类型）
#   gwat       :: 来自上部的净水输入（mm H₂O/s，T 类型）
# =============================================================================
function surface_runoff_simtop(N::Int,
  Ksat::Vector{T}, F_max::T, f::T,
  ice_frac::Vector{T}, zwt::T, gwat::T) where T<:Real

  # 饱和面积比例：当水位较浅时，fsat 越接近 fsatmax；当水位增深时按指数衰减
  F_sat = F_max * min(1.0, exp(-T(2.0) * f * zwt))

  n = min(3, N) # 仅考虑前3层
  I_candidates = [T(10.0)^(-T(6.0) * ice_frac[i]) * Ksat[i] for i in 1:n]
  I_max = minimum(I_candidates)

  ## 若第一层有效孔隙率低于阈值 wimp，则入渗能力置为 0
  # eff_porosity[1] < wimp && (I_max = 0.0) # effective porosity = porosity - vol_ice
  rsur_se = F_sat * max(0.0, gwat)                  # 饱和区, Niu 2007, JGR-A, Eq.1
  rsur_ie = (1.0 - F_sat) * max(0.0, gwat - I_max)  # 下渗

  rsur = rsur_se + rsur_ie
  return rsur, rsur_se, rsur_ie
end

# =============================================================================
# 函数：subsurface_runoff_simtop
#
# 功能：计算地下径流（单位：mm H₂O/s），考虑冰冻对渗流的阻碍作用。
#
# 参数说明：
#   N    :: 土壤层数
#   icefrac    :: 各层冰的比例，向量
#   dz_soisno  :: 各层厚度（m），向量
#   zi_soisno  :: 各层界面位置（m），向量（长度应为 N+1，对应 Fortran 中的 0:N）
#   zwt        :: 水位深度（m，T 类型）
# =============================================================================
function subsurface_runoff_simtop(N::Int, icefrac::Vector{T},
  dz_soisno::Vector{T}, zi_soisno::Vector{T},
  zwt::T) where T<:Real

  dzmm = dz_soisno .* T(1000.0) # [m] to [mm]

  # 确定水位上方土层的索引（变量 jwt）
  # 如果水位位于最上层，则 jwt 设为 0；否则找到第一层使得 zwt 小于界面位置
  jwt = N
  for j in 1:N
    if zwt <= zi_soisno[j]
      jwt = j - 1
      break
    end
  end
  jwt = find_jwt(zi_soisno, zwt)

  # 累加从确定层开始的各层厚度和冰分数加权值
  dzsum = 0.0
  icefracsum = 0.0

  start_index = max(jwt, 1) # 如果 jwt 为 0，则从第 1 层开始累计
  for j in start_index:N
    dzsum += dzmm[j]
    icefracsum += icefrac[j] * dzmm[j]
  end

  # 计算冰冻抑制因子：
  # fracice_rsub = [exp(-3*(1-冰加权分数))-exp(-3)] / [1-exp(-3)]
  frac_ice = max(0.0, exp(-T(3.0) * (1.0 - (icefracsum / dzsum))) - exp(-T(3.0))) /
             (1.0 - exp(-T(3.0)))
  imped = max(0.0, 1.0 - frac_ice) # 阻碍因子
  Rsb_max = imped * T(5.5e-3)
  rsubst = Rsb_max * exp(-T(2.5) * zwt)
  return rsubst
end

# =============================================================================
# 函数：runoff_xinanjiang
#
# 功能：采用新安江方法计算径流，返回一个二元组 (rsur, rsubst)，其中：
#       - rsur    ：地表径流（单位：mm H₂O/s）
#       - rsubst ：地下径流（此方法中固定为 0）
#
# 参数说明：
#   N      :: 土壤层数
#   dz_soisno    :: 各层厚度（m），向量
#   eff_porosity :: 各层有效孔隙率，向量
#   vol_liq      :: 各层液态水体积分数，向量
#   topostd      :: 地形高程标准差（m，T 类型）
#   gwat         :: 上部水输入（mm H₂O/s，T 类型）
#   deltim       :: 时间步长（s，T 类型）
# =============================================================================
function runoff_xinanjiang(N::Int, dz_soisno::Vector{T},
  eff_porosity::Vector{T}, vol_liq::Vector{T},
  topostd::T, gwat::T, deltim::T) where T<:Real
  sigmin = T(100.0)
  sigmax = T(1000.0)

  # 将水输入量转换为以 m 为单位（watin）
  watin = gwat * deltim / T(1000.0)

  if watin <= 0.0
    rsur = 0.0
    rsubst = 0.0
    return rsur, rsubst
  else
    # 计算地形因子 btopo，并限制其范围在 [0.01, 0.5] 内
    btopo = (topostd - sigmin) / (topostd + sigmax)
    btopo = min(max(btopo, T(0.01)), T(0.5))

    # 对前 6 层土壤累计（如果土层数小于 6，则取 N）
    n_layer = min(6, N)
    w_int = sum(vol_liq[1:n_layer] .* dz_soisno[1:n_layer])
    wsat_int = sum(eff_porosity[1:n_layer] .* dz_soisno[1:n_layer])

    # 计算中间变量 wtmp，并由此得到入渗量 infil
    wtmp = (1.0 - w_int / wsat_int)^(1.0 / (btopo + 1.0)) - watin / ((btopo + 1.0) * wsat_int)
    infil = wsat_int - w_int - wsat_int * (max(0.0, wtmp))^(btopo + 1.0)
    infil = min(infil, watin)

    # 剩余水量作为地表径流 rsur（转换回 mm H₂O/s）
    rsur = (watin - infil) * T(1000.0) / deltim
    rsubst = 0.0
    return rsur, rsubst
  end
end

# =============================================================================
# 函数：runoff_simplevic
#
# 功能：采用简化版 VIC 模型计算径流，返回二元组 (rsur, rsubst)，其中：
#       - rsur    ：地表径流（单位：mm H₂O/s）
#       - rsubst ：地下径流（此方法中固定为 0）
#
# 参数说明：
#   N      :: 土壤层数
#   dz_soisno    :: 各层厚度（m），向量
#   eff_porosity :: 各层有效孔隙率，向量
#   vol_liq      :: 各层液态水体积分数，向量
#   BVIC         :: VIC 模型参数（T 类型）
#   gwat         :: 上部水输入（mm H₂O/s，T 类型）
#   deltim       :: 时间步长（s，T 类型）
# =============================================================================
function runoff_simplevic(N::Int, dz_soisno::Vector{T},
  eff_porosity::Vector{T}, vol_liq::Vector{T},
  BVIC::T, gwat::T, dt::T) where T<:Real
  # 将水输入量转换为以 m 为单位
  watin = gwat * dt / T(1000.0)

  if watin <= 0.0
    rsur = 0.0
    rsubst = 0.0
    return rsur, rsubst
  else
    # 对前 6 层土壤累计（取 N 与 6 中的较小值）
    n_layer = min(6, N)
    w_int = sum(vol_liq[1:n_layer] .* dz_soisno[1:n_layer])
    wsat_int = sum(eff_porosity[1:n_layer] .* dz_soisno[1:n_layer])

    # 计算土壤饱和分数（SoilSaturateFrac）
    InfilExpFac = BVIC / (1.0 + BVIC)
    SoilSaturateFrac = 1.0 - (max(0.0, 1.0 - (w_int / wsat_int)))^InfilExpFac
    SoilSaturateFrac = max(0.0, min(1.0, SoilSaturateFrac))

    # 计算最大水深和初始水深
    WaterDepthMax = (1.0 + BVIC) * wsat_int
    WaterDepthInit = WaterDepthMax * (1.0 - (1.0 - SoilSaturateFrac)^(1.0 / BVIC))

    # 根据条件判断计算地表径流
    if WaterDepthMax <= 0.0
      Rs = watin
    elseif (WaterDepthInit + watin) > WaterDepthMax
      Rs = watin - wsat_int + w_int
    else
      InfilVarTmp = 1.0 - ((WaterDepthInit + watin) / WaterDepthMax)
      Rs = watin - wsat_int + w_int + wsat_int * (InfilVarTmp^(1.0 + BVIC))
    end

    # 限制径流值在 [0, watin] 范围内
    Rs = max(0.0, min(Rs, watin))
    # 计算入渗量（辅助计算，不单独输出）
    infil = watin - Rs
    rsur = Rs * T(1000.0) / dt
    rsubst = 0.0
    return rsur, rsubst
  end
end

end # module Runoff
