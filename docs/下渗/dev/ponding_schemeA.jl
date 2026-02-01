"""
# 方案A：基于蓄水容量的地表积水处理方案

## 物理机制
基于蓄满产流机制，定义最大地表蓄水深度 `h_pond_max`，代表微地形蓄水容量。
当积水深度超过此阈值时，超出部分立即产生径流。

## 适用场景
- 平坦地区
- 农田
- 城市下垫面

## 优点
1. 物理意义明确：h_pond_max代表微地形蓄水容量
2. 与Richards方程耦合：h₀直接影响下渗能力
3. 符合蓄满产流机制
4. 计算简单高效

## 参数
- `h_pond_max`: 最大地表蓄水深度 [cm]，典型值 0.1-2.0 cm
  - 平坦农田：0.3-0.5 cm
  - 城市不透水面：0.1-0.3 cm
  - 有一定起伏的地形：0.5-1.0 cm

## 参考文献
见 docs/下渗/Infiltration_final.typ
"""


"""
    ponding_schemeA(h_pond_current, P, inf, dt; h_pond_max=0.5)

方案A：基于蓄水容量的地表积水处理

# Arguments
- `h_pond_current`: 当前积水深度 [cm]
- `P`: 降水强度 [cm/h]
- `inf`: 实际下渗速率 [cm/h]（正值表示下渗）
- `dt`: 时间步长 [h]
- `h_pond_max`: 最大蓄水深度 [cm]，默认 0.5 cm

# Returns
- `runoff`: 产生的径流深度 [cm]
- `h_pond_new`: 更新后的积水深度 [cm]

# Example
```julia
# 初始积水 0.2 cm，降雨 5 cm/h，下渗 2 cm/h，时间步长 0.1 h
runoff, h_pond = ponding_schemeA(0.2, 5.0, 2.0, 0.1; h_pond_max=0.5)
# 积水增量 = (5.0 - 2.0) * 0.1 = 0.3 cm
# 新积水 = 0.2 + 0.3 = 0.5 cm (未超过阈值)
# 径流 = 0.0 cm
```
"""
function ponding_schemeA(h_pond_current::T, P::T, inf::T, dt::T;
                        h_pond_max::T=T(0.5), ignored...) where {T<:Real}
  # 地表水量平衡
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


"""
    ponding_schemeA!(soil, h_pond_current, P, inf; h_pond_max=0.5)

方案A的原位修改版本（修改soil对象）

# Arguments
- `soil`: Soil对象
- `h_pond_current`: 当前积水深度 [cm]
- `P`: 降水强度 [cm/h]
- `inf`: 实际下渗速率 [cm/h]（正值表示下渗）
- `h_pond_max`: 最大蓄水深度 [cm]，默认 0.5 cm

# Returns
- `runoff`: 产生的径流深度 [cm]
- `h_pond_new`: 更新后的积水深度 [cm]

# Notes
此函数使用 `soil.dt` 作为时间步长（单位：秒），自动转换为小时。
"""
function ponding_schemeA!(soil, h_pond_current::T, P::T, inf::T;
                         h_pond_max::T=T(0.5)) where {T<:Real}
  dt_h = soil.dt / 3600.0  # 转换为小时
  return ponding_schemeA(h_pond_current, P, inf, dt_h; h_pond_max=h_pond_max)
end
