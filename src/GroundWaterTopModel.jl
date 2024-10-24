# module GroundWaterTopModelMod
# 模拟地下水流动和基于TOPMODEL的地下径流 (Niu et al., 2007)

mutable struct NoahmpType{FT}
  n::Int                      # 土壤层数
  dt::FT                      # 土壤时间步长 [s]
  z::Vector{FT}               # 土壤层底部深度 [m], 向下为负
  F_frz_max::FT                  # 最大土壤不透水部分 [-]
  θ_ice::Vector{FT}                    # 土壤冰含量 [m3/m3]
  K::Vector{FT}        # 土壤水力传导率 [m/s]
  θ_sat::Vector{FT}            # 土壤饱和水分含量 [m3/m3]
  λ::FT                      # 网格平均地形指数 [-]
  Ψ_sat::Vector{FT}        # 饱和基质势 [m]
  SoilExpCoeffB::Vector{FT}              # 土壤B参数 [-]
  SpecYieldGw::FT                        # 特定产水量 [-]
  MicroPoreContent::FT                   # 微孔含量 (0.0-1.0)
  Ksat::Vector{FT}     # 饱和土壤水力传导率 [m/s]
  θ_liq::Vector{FT}               # 土壤液态水含量 [m3/m3]
  WaterTableDepth::FT                    # 地下水位深度 [m]
  WaterStorageAquifer::FT                # 含水层水储量 [mm]
  WaterStorageSoilAqf::FT                # 含水层和饱和土壤的水储量 [mm]
  RunoffDecayFac::FT                     # 径流衰减因子 (1/m)
  BaseflowCoeff::FT                      # 底流系数 [mm/s]
  RechargeGw::FT                         # 地下水补给速率 [mm/s]
  DischargeGw::FT                        # 地下水排放速率 [mm/s]
end

function GroundWaterTopModel!(noahmp::NoahmpType{FT}) where {FT}
  (; n, dt, z, F_frz_max, θ_ice, K,
    θ_sat, λ, Ψ_sat, SoilExpCoeffB, Sy, MicroPoreContent,
    Ksat, θ_liq, hgw, Wa, S_SoilAqf,
    f, R_sb_max, Q, R_sb) = noahmp

  dz_mm = zeros(FT, n)    # 每层土壤的厚度
  h₊ₕ_m = zeros(FT, n)   # 每层的中间位置
  _θ_liq = zeros(FT, n)
  porosity = zeros(FT, n)
  _K = zeros(FT, n)
  θ = zeros(FT, n)                 # SoilMoisture

  # 计算土壤层厚度和节点深度
  dz_mm[1] = -z[1] * 1.0e3
  for i in 2:n
    dz_mm[i] = 1.0e3 * (z[i-1] - z[i])
  end

  h₊ₕ_m[1] = -z[1] / 2.0
  for i in 2:n
    h₊ₕ_m[i] = -z[i-1] + 0.5 * (z[i-1] - z[i])
  end

  # 将体积土壤水分转换为质量
  for i in 1:n
    θ[i] = θ_liq[i] + θ_ice[i]
    _θ_liq[i] = θ_liq[i] * dz_mm[i]
    porosity[i] = max(0.01, θ_sat[i] - θ_ice[i]) # effective porosity
    _K[i] = 1.0e3 * K[i]
  end

  # 查找第一层未饱和层的索引
  j = n # the unsaturated layer, above the water table
  for i in 2:n
    if hgw <= -z[i]
      j = i - 1
      break
    end
  end

  # 计算地下水排放
  f = SoilExpCoeffB[j] / 3.0 # runoff decay factor
  R_sb_max = _K[j] * 1.0e3 * exp(3.0)
  R_sb = (1.0 - F_frz_max) * R_sb_max * exp(-λ - f * hgw)

  # 计算未饱和层的基质势
  perc_unsat = clamp(θ[j] / θ_sat[j], 0.01, 1.0)

  Ψ_Frz = -Ψ_sat[j] * 1000.0 * perc_unsat^(-SoilExpCoeffB[j])
  Ψ_Frz = max(-120000.0, MicroPoreContent * Ψ_Frz)

  # 计算地下水补给速率
  Ka = 2.0 * (_K[j] * Ksat[j] * 1.0e3) /
       (_K[j] + Ksat[j] * 1.0e3)
  zgw = -hgw * 1.0e3
  z_head = Ψ_Frz - h₊ₕ_m[j] * 1.0e3
  Q = -Ka * (zgw - z_head) /
      ((hgw - h₊ₕ_m[j]) * 1.0e3)
  Q = max(-10.0 / dt, min(10.0 / dt, Q))

  # 更新含水层储水量和地下水位
  S_SoilAqf += (Q - R_sb) * dt

  if j == n
    Wa += (Q - R_sb) * dt
    S_SoilAqf = Wa
    hgw = (-z[n] + 25.0) - Wa / 1000.0 / Sy
    _θ_liq[n] -= Q * dt
    _θ_liq[n] += max(0.0, Wa - 5000.0)
    Wa = min(Wa, 5000.0)
  else
    # 更新其他层的水位
    WaterFillPore = 0.0
    for i in (j+2):n
      WaterFillPore += porosity[i] * dz_mm[i]
    end
    hgw = -z[j+1] - (S_SoilAqf - Sy * 1000.0 * 25.0 - WaterFillPore) / (porosity[j+1]) / 1000.0
  end

  WatConductAcc = 0.0
  for i in 1:n
    WatConductAcc += _K[i] * dz_mm[i]
  end
  for i in 1:n
    _θ_liq[i] -= R_sb * dt * _K[i] * dz_mm[i] / WatConductAcc
  end

  hgw = max(1.5, hgw)

  # 限制SoilLiqTmp不小于SoilMoistureMin
  θ_min = 0.01
  for i in 1:n-1
    if _θ_liq[i] < 0.0
      θ_excess = θ_min - _θ_liq[i]
    else
      θ_excess = 0.0
    end
    _θ_liq[i] += θ_excess
    _θ_liq[i+1] -= θ_excess
  end

  if _θ_liq[n] < θ_min
    θ_excess = θ_min - _θ_liq[n]
  else
    θ_excess = 0.0
  end
  _θ_liq[n] += θ_excess
  Wa -= θ_excess
  S_SoilAqf -= θ_excess

  # 更新土壤水分
  for i in 1:n
    θ_liq[i] = _θ_liq[i] / dz_mm[i]
  end
end
