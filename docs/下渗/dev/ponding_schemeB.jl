"""
# 方案B：基于坡面汇流的地表积水处理方案

## 物理机制
基于曼宁公式（Manning's equation）计算坡面径流速度，考虑地形坡度、
地表糙度和积水深度的非线性影响。

## 理论基础

### 曼宁公式
```
v = (1/n) * h^(2/3) * S₀^(1/2)
```

### 单位宽度流量
```
q = v * h = (1/n) * h^(5/3) * S₀^(1/2)
```

### 排水速率
```
drain_rate = q/L = (1/(n*L)) * h^(5/3) * S₀^(1/2)
```

其中 L 是坡面特征长度，来源于连续性方程的简化：
```
∂h/∂t = -∂q/∂x ≈ -q/L
```

## 适用场景
- 坡地
- 山区
- 自然流域

## 优点
1. 考虑地形坡度和地表糙度
2. 积水深度越大，排水越快（h^(5/3) 非线性关系）
3. 物理过程清晰

## 参数

### 曼宁粗糙系数 n [s/m^(1/3)]
- 裸地/光滑表面：0.01-0.05
- 农田：0.10-0.20
- 草地：0.15-0.40
- 森林/灌木：0.30-0.80

### 坡度 S₀ [-]
- 平坦：0.001-0.01 (0.1%-1%)
- 缓坡：0.01-0.05 (1%-5%)
- 中坡：0.05-0.15 (5%-15%)
- 陡坡：0.15-0.30 (15%-30%)

### 特征长度 L [m]
单站点模拟推荐值：
- 农田/草地：30-50 m
- 森林：50-100 m
- 城市：10-30 m
- 快速测试：50 m

### 起始产流深度 h_pond_min [cm]
代表微地形蓄水，典型值 0.1-0.5 cm

## 参考文献
见 docs/下渗/Infiltration_final.typ
"""


"""
    ponding_schemeB(h_pond_current, P, inf, dt;
                    S0=0.05, n_manning=0.15, L_char=50.0, h_pond_min=0.1)

方案B：基于坡面汇流的地表积水处理（曼宁公式）

# Arguments
- `h_pond_current`: 当前积水深度 [cm]
- `P`: 降水强度 [cm/h]
- `inf`: 实际下渗速率 [cm/h]（正值表示下渗）
- `dt`: 时间步长 [h]
- `S0`: 地面坡度 [-]，默认 0.05 (5%)
- `n_manning`: 曼宁粗糙系数 [s/m^(1/3)]，默认 0.15
- `L_char`: 特征汇流长度 [m]，默认 50.0
- `h_pond_min`: 起始产流深度 [cm]，默认 0.1

# Returns
- `runoff`: 产生的径流深度 [cm]
- `h_pond_new`: 更新后的积水深度 [cm]

# Physics
基于曼宁公式：
- 流速：v = (1/n) * h^(2/3) * S₀^(1/2) [m/s]
- 单位流量：q = v * h = (1/n) * h^(5/3) * S₀^(1/2) [m²/s]
- 排水速率：drain_rate = q/L [m/s] → [cm/h]

# Example
```julia
# 坡度5%，农田糙度，坡长50m
runoff, h_pond = ponding_schemeB(1.0, 5.0, 2.0, 0.1;
                                 S0=0.05, n_manning=0.15, L_char=50.0)
```

# Notes
- h_pond_min 以下不产流（微地形蓄水）
- 排水速率随 h^(5/3) 增长（非线性强）
- 单位转换：内部计算使用 [m]，输入输出使用 [cm]
"""
function ponding_schemeB(h_pond_current::T, P::T, inf::T, dt::T;
                        S0::T=T(0.05),
                        n_manning::T=T(0.15),
                        L_char::T=T(50.0),
                        h_pond_min::T=T(0.1), ignored...) where {T<:Real}
  # 地表水量平衡
  h_pond = h_pond_current + (P - inf) * dt  # [cm]

  # 产流计算（简化运动波）
  if h_pond > h_pond_min
    # 计算超过起始产流深度的部分
    h_excess = h_pond - h_pond_min  # [cm]
    h_m = h_excess / 100.0  # 转换为 [m]

    # 曼宁公式：v = (1/n) * h^(2/3) * S^(1/2)
    # 假设宽浅水流：水力半径 R ≈ h
    v = (1.0 / n_manning) * (h_m^(2/3)) * sqrt(S0)  # [m/s]

    # 单位宽度流量 q = v * h
    q = v * h_m  # [m²/s]
    q_cm_h = q * 3600.0 * 100.0  # [m²/s] → [cm²/h]

    # 排水速率：drain_rate = q / L
    drain_rate = q_cm_h / L_char  # [cm/h]

    # 更新积水深度（不能降到 h_pond_min 以下）
    drainage = min(drain_rate * dt, h_excess)
    runoff = drainage
    h_pond_new = h_pond - drainage
  else
    runoff = 0.0
    h_pond_new = max(h_pond, 0.0)
  end

  return runoff, h_pond_new
end


"""
    ponding_schemeB!(soil, h_pond_current, P, inf;
                     S0=0.05, n_manning=0.15, L_char=50.0, h_pond_min=0.1)

方案B的原位修改版本（修改soil对象）

# Arguments
- `soil`: Soil对象
- `h_pond_current`: 当前积水深度 [cm]
- `P`: 降水强度 [cm/h]
- `inf`: 实际下渗速率 [cm/h]（正值表示下渗）
- `S0`: 地面坡度 [-]，默认 0.05
- `n_manning`: 曼宁粗糙系数 [s/m^(1/3)]，默认 0.15
- `L_char`: 特征汇流长度 [m]，默认 50.0
- `h_pond_min`: 起始产流深度 [cm]，默认 0.1

# Returns
- `runoff`: 产生的径流深度 [cm]
- `h_pond_new`: 更新后的积水深度 [cm]

# Notes
此函数使用 `soil.dt` 作为时间步长（单位：秒），自动转换为小时。
"""
function ponding_schemeB!(soil, h_pond_current::T, P::T, inf::T;
                         S0::T=T(0.05),
                         n_manning::T=T(0.15),
                         L_char::T=T(50.0),
                         h_pond_min::T=T(0.1)) where {T<:Real}
  dt_h = soil.dt / 3600.0  # 转换为小时
  return ponding_schemeB(h_pond_current, P, inf, dt_h;
                        S0=S0, n_manning=n_manning, L_char=L_char, h_pond_min=h_pond_min)
end
