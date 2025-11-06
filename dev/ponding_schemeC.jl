"""
# 方案C：混合方案（推荐）

## 物理机制
综合考虑微地形蓄水和坡面汇流两个过程，采用三阶段机制：

1. **微地形蓄水阶段** (h < h_pond_min)：
   - 降水填充微地形凹陷，无径流产生

2. **坡面汇流阶段** (h_pond_min < h < h_pond_max)：
   - 采用曼宁公式计算径流速率
   - 径流速率随积水深度非线性增加（h^(5/3)）

3. **超蓄阶段** (h > h_pond_max)：
   - 超过地表最大调蓄能力
   - 超出部分立即产流

## 适用场景
通用场景，特别适合：
- 需要综合考虑地形和糙度的模拟
- 有一定坡度但不太陡的地区
- 参数可通过观测率定的场景

## 优点
1. 物理机制完整，分阶段处理不同情况
2. 综合方案A和B的优点
3. 参数物理意义明确，可通过观测率定
4. 适应性强

## 参数
- `h_pond_min`: 微地形蓄水容量 [cm]，典型值 0.1-0.5 cm
- `h_pond_max`: 最大允许积水深度 [cm]，典型值 1.0-5.0 cm
- `S0`: 坡度 [-]，典型值 0.001-0.10
- `n_manning`: 曼宁系数 [s/m^(1/3)]，见方案B
- `L_char`: 汇流特征长度 [m]，典型值 30-100 m

## 参考文献
见 docs/下渗/Infiltration_final.typ
"""


"""
    ponding_schemeC(h_pond_current, P, inf, dt;
                    h_pond_min=0.1, h_pond_max=2.0,
                    S0=0.01, n_manning=0.15, L_char=100.0)

方案C：混合方案（微地形蓄水 + 坡面汇流 + 超蓄产流）

# Arguments
- `h_pond_current`: 当前积水深度 [cm]
- `P`: 降水强度 [cm/h]
- `inf`: 实际下渗速率 [cm/h]（正值表示下渗）
- `dt`: 时间步长 [h]
- `h_pond_min`: 微地形蓄水容量 [cm]，默认 0.1
- `h_pond_max`: 最大允许积水深度 [cm]，默认 2.0
- `S0`: 坡度 [-]，默认 0.01 (1%)
- `n_manning`: 曼宁系数 [s/m^(1/3)]，默认 0.15
- `L_char`: 汇流特征长度 [m]，默认 100.0

# Returns
- `runoff`: 产生的径流深度 [cm]
- `h_pond_new`: 更新后的积水深度 [cm]

# Physics

## 三阶段机制

### 阶段1：h ≤ h_pond_min
微地形蓄水，无径流：
```
runoff = 0
```

### 阶段2：h_pond_min < h ≤ h_pond_max
坡面汇流（曼宁公式）：
```
drain_rate = (1/(n*L)) * (h - h_pond_min)^(5/3) * S₀^(1/2)
runoff = drain_rate * dt
```

### 阶段3：h > h_pond_max
超蓄产流（立即排出）：
```
runoff = h - h_pond_max
```

# Example
```julia
# 微坡度，中等糙度，较长坡面
runoff, h_pond = ponding_schemeC(0.5, 10.0, 3.0, 0.1;
                                 h_pond_min=0.1, h_pond_max=2.0,
                                 S0=0.01, n_manning=0.15, L_char=100.0)
```

# Notes
- 结合了方案A的蓄水容量和方案B的曼宁公式
- 三阶段平滑过渡
- 对于平坦地区（S0很小），主要依靠h_pond_max控制
- 对于坡地（S0较大），坡面汇流起主导作用
"""
function ponding_schemeC(h_pond_current::T, P::T, inf::T, dt::T;
                        h_pond_min::T=T(0.1),
                        h_pond_max::T=T(2.0),
                        S0::T=T(0.01),
                        n_manning::T=T(0.15),
                        L_char::T=T(100.0), ignored...) where {T<:Real}
  # 地表水量平衡
  h_pond = h_pond_current + (P - inf) * dt  # [cm]

  # 分段处理
  if h_pond <= h_pond_min
    # 阶段1：微地形蓄水阶段，无径流
    runoff = 0.0
    h_pond_new = max(h_pond, 0.0)

  elseif h_pond <= h_pond_max
    # 阶段2：坡面汇流阶段，运动波方程
    h_excess = h_pond - h_pond_min
    h_m = h_excess / 100.0  # [cm] → [m]

    # 曼宁公式
    v = (1.0 / n_manning) * (h_m^(2/3)) * sqrt(max(S0, 1e-6))  # [m/s]
    q = v * h_m  # [m²/s]
    q_cm_h = q * 3600.0 * 100.0  # [cm²/h]
    drain_rate = q_cm_h / L_char  # [cm/h]

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


"""
    ponding_schemeC!(soil, h_pond_current, P, inf;
                     h_pond_min=0.1, h_pond_max=2.0,
                     S0=0.01, n_manning=0.15, L_char=100.0)

方案C的原位修改版本（修改soil对象）

# Arguments
- `soil`: Soil对象
- `h_pond_current`: 当前积水深度 [cm]
- `P`: 降水强度 [cm/h]
- `inf`: 实际下渗速率 [cm/h]（正值表示下渗）
- `h_pond_min`: 微地形蓄水容量 [cm]，默认 0.1
- `h_pond_max`: 最大允许积水深度 [cm]，默认 2.0
- `S0`: 坡度 [-]，默认 0.01
- `n_manning`: 曼宁系数 [s/m^(1/3)]，默认 0.15
- `L_char`: 汇流特征长度 [m]，默认 100.0

# Returns
- `runoff`: 产生的径流深度 [cm]
- `h_pond_new`: 更新后的积水深度 [cm]

# Notes
此函数使用 `soil.dt` 作为时间步长（单位：秒），自动转换为小时。
"""
function ponding_schemeC!(soil, h_pond_current::T, P::T, inf::T;
                         h_pond_min::T=T(0.1),
                         h_pond_max::T=T(2.0),
                         S0::T=T(0.01),
                         n_manning::T=T(0.15),
                         L_char::T=T(100.0)) where {T<:Real}
  dt_h = soil.dt / 3600.0  # 转换为小时
  return ponding_schemeC(h_pond_current, P, inf, dt_h;
                        h_pond_min=h_pond_min, h_pond_max=h_pond_max,
                        S0=S0, n_manning=n_manning, L_char=L_char)
end
