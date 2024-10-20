# module GroundWaterTopModelMod
# 模拟地下水流动和基于TOPMODEL的地下径流 (Niu et al., 2007)

mutable struct NoahmpType{FT}
  n::Int                      # 土壤层数
  dt::FT                      # 土壤时间步长 [s]
  z::Vector{FT}               # 土壤层底部深度 [m], 向下为负
  SoilImpervFracMax::FT                  # 最大土壤不透水部分 [-]
  SoilIce::Vector{FT}                    # 土壤冰含量 [m3/m3]
  SoilWatConductivity::Vector{FT}        # 土壤水力传导率 [m/s]
  Θ_sat::Vector{FT}            # 土壤饱和水分含量 [m3/m3]
  GridTopoIndex::FT                      # 网格平均地形指数 [-]
  SoilMatPotentialSat::Vector{FT}        # 饱和基质势 [m]
  SoilExpCoeffB::Vector{FT}              # 土壤B参数 [-]
  SpecYieldGw::FT                        # 特定产水量 [-]
  MicroPoreContent::FT                   # 微孔含量 (0.0-1.0)
  SoilWatConductivitySat::Vector{FT}     # 饱和土壤水力传导率 [m/s]
  SoilLiqWater::Vector{FT}               # 土壤液态水含量 [m3/m3]
  WaterTableDepth::FT                    # 地下水位深度 [m]
  WaterStorageAquifer::FT                # 含水层水储量 [mm]
  WaterStorageSoilAqf::FT                # 含水层和饱和土壤的水储量 [mm]
  RunoffDecayFac::FT                     # 径流衰减因子 (1/m)
  BaseflowCoeff::FT                      # 底流系数 [mm/s]
  RechargeGw::FT                         # 地下水补给速率 [mm/s]
  DischargeGw::FT                        # 地下水排放速率 [mm/s]
end

function GroundWaterTopModel!(noahmp::NoahmpType{FT}) where {FT}
  (; n, dt, z, SoilImpervFracMax, SoilIce, SoilWatConductivity,
   Θ_sat, GridTopoIndex, SoilMatPotentialSat, SoilExpCoeffB, SpecYieldGw, MicroPoreContent,
   SoilWatConductivitySat, SoilLiqWater, WaterTableDepth, WaterStorageAquifer, WaterStorageSoilAqf,
   RunoffDecayFac, BaseflowCoeff, RechargeGw, DischargeGw) = noahmp
   
  # 初始化
  dz = zeros(FT, n)    # 每层土壤的厚度
  h₊ₕ = zeros(FT, n)   # 每层的中间位置
  
  SoilLiqTmp = zeros(FT, n)
  SoilEffPorosity = zeros(FT, n)
  SoilWatConductTmp = zeros(FT, n)
  Θ = zeros(FT, n)

  # 计算土壤层厚度和节点深度
  dz[1] = -z[1] * 1.0e3
  for i in 2:n
    dz[i] = 1.0e3 * (z[i-1] - z[i])
  end

  h₊ₕ[1] = -z[1] / 2.0
  for i in 2:n
    h₊ₕ[i] = -z[i-1] + 0.5 * (z[i-1] - z[i])
  end

  # 将体积土壤水分转换为质量
  for i in 1:n
    Θ[i] = SoilLiqWater[i] + SoilIce[i]
    SoilLiqTmp[i] = SoilLiqWater[i] * dz[i]
    SoilEffPorosity[i] = max(0.01, Θ_sat[i] - SoilIce[i])
    SoilWatConductTmp[i] = 1.0e3 * SoilWatConductivity[i]
  end

  # 查找第一层未饱和层的索引
  IndUnsatSoil = n
  for i in 2:n
    if WaterTableDepth <= -z[i]
      IndUnsatSoil = i - 1
      break
    end
  end

  # 计算地下水排放
  RunoffDecayFac = SoilExpCoeffB[IndUnsatSoil] / 3.0
  BaseflowCoeff = SoilWatConductTmp[IndUnsatSoil] * 1.0e3 * exp(3.0)
  DischargeGw = (1.0 - SoilImpervFracMax) * BaseflowCoeff * exp(-GridTopoIndex) *
                exp(-RunoffDecayFac * WaterTableDepth)

  # 计算未饱和层的基质势
  SatDegUnsatSoil = min(1.0, Θ[IndUnsatSoil] / Θ_sat[IndUnsatSoil])
  SatDegUnsatSoil = max(SatDegUnsatSoil, 0.01)
  SoilMatPotFrz = -SoilMatPotentialSat[IndUnsatSoil] * 1000.0 * SatDegUnsatSoil^(-SoilExpCoeffB[IndUnsatSoil])
  SoilMatPotFrz = max(-120000.0, MicroPoreContent * SoilMatPotFrz)

  # 计算地下水补给速率
  AquiferWatConduct = 2.0 * (SoilWatConductTmp[IndUnsatSoil] * SoilWatConductivitySat[IndUnsatSoil] * 1.0e3) /
                      (SoilWatConductTmp[IndUnsatSoil] + SoilWatConductivitySat[IndUnsatSoil] * 1.0e3)
  WaterHeadTbl = -WaterTableDepth * 1.0e3
  WaterHead = SoilMatPotFrz - h₊ₕ[IndUnsatSoil] * 1.0e3
  RechargeGw = -AquiferWatConduct * (WaterHeadTbl - WaterHead) /
               ((WaterTableDepth - h₊ₕ[IndUnsatSoil]) * 1.0e3)
  RechargeGw = max(-10.0 / dt, min(10.0 / dt, RechargeGw))

  # 更新含水层储水量和地下水位
  WaterStorageSoilAqf += (RechargeGw - DischargeGw) * dt
  if IndUnsatSoil == n
    WaterStorageAquifer += (RechargeGw - DischargeGw) * dt
    WaterStorageSoilAqf = WaterStorageAquifer
    WaterTableDepth = (-z[n] + 25.0) - WaterStorageAquifer / 1000.0 / SpecYieldGw
    SoilLiqTmp[n] -= RechargeGw * dt
    SoilLiqTmp[n] += max(0.0, WaterStorageAquifer - 5000.0)
    WaterStorageAquifer = min(WaterStorageAquifer, 5000.0)
  else
    # 更新其他层的水位
    WaterFillPore = 0.0
    for i in (IndUnsatSoil+2):n
      WaterFillPore += SoilEffPorosity[i] * dz[i]
    end
    WaterTableDepth = -z[IndUnsatSoil+1] - (WaterStorageSoilAqf - SpecYieldGw * 1000.0 * 25.0 -
                                                         WaterFillPore) / (SoilEffPorosity[IndUnsatSoil+1]) / 1000.0
  end

  WatConductAcc = 0.0
  for i in 1:n
    WatConductAcc += SoilWatConductTmp[i] * dz[i]
  end
  for i in 1:n
    SoilLiqTmp[i] -= DischargeGw * dt * SoilWatConductTmp[i] * dz[i] / WatConductAcc
  end

  WaterTableDepth = max(1.5, WaterTableDepth)

  # 限制SoilLiqTmp不小于SoilMoistureMin
  Θ_min = 0.01
  for i in 1:n-1
    if SoilLiqTmp[i] < 0.0
      WaterExcessSat = Θ_min - SoilLiqTmp[i]
    else
      WaterExcessSat = 0.0
    end
    SoilLiqTmp[i] += WaterExcessSat
    SoilLiqTmp[i+1] -= WaterExcessSat
  end

  if SoilLiqTmp[n] < Θ_min
    WaterExcessSat = Θ_min - SoilLiqTmp[n]
  else
    WaterExcessSat = 0.0
  end
  SoilLiqTmp[n] += WaterExcessSat
  WaterStorageAquifer -= WaterExcessSat
  WaterStorageSoilAqf -= WaterExcessSat

  # 更新土壤水分
  for i in 1:n
    SoilLiqWater[i] = SoilLiqTmp[i] / dz[i]
  end
end
