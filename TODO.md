# 下渗模块开发计划

## 方案 1: 降雨-下渗耦合模块（推荐）

在 `Equation_Richards.jl` 基础上添加自适应边界条件：

```julia
function infiltration_boundary!(soil::Soil, P::Float64, dt::Float64)
  """
  P: precipitation rate [cm/h]
  自动切换边界条件：
  - 无积水时: Q0 = -P (通量边界)
  - 积水时: ψ0 = 0 (饱和边界)
  """

  # 计算潜在下渗能力
  K_surf = soil.K[soil.ibeg]
  ψ_surf = soil.ψ[soil.ibeg]

  # Case 1: 降雨强度 < 下渗能力 → 通量边界
  if P <= K_surf && soil.uex <= 0
    method = "Q0"
    soil.Q0 = -P
    soil.ψ0 = NaN

  # Case 2: 降雨强度 > 下渗能力 → 可能积水
  else
    # 尝试通量边界
    Q0_potential = cal_Q!(soil; Q0=-P, method="Q0")

    # 检查表层是否达到饱和
    θ_new = soil.θ[ibeg] + (Q0_potential - soil.Q[ibeg]) * dt/3600 / soil.Δz[ibeg]

    if θ_new >= soil.param.θ_sat[ibeg]
      # 切换到饱和边界
      method = "ψ0"
      soil.ψ0 = 0.0  # 饱和
      soil.uex += (P - actual_infiltration) * dt/3600  # 积水
    else
      method = "Q0"
      soil.Q0 = -P
    end
  end

  return method
end
```

## 方案 2: Green-Ampt 简化模型

作为快速替代方案（计算效率高）:

```julia
mutable struct GreenAmpt{T}
  Ksat::T          # 饱和导水率 [cm/h]
  θi::T            # 初始含水量
  θs::T            # 饱和含水量
  ψf::T            # 湿润锋吸力 [cm]
  F::T = 0.0       # 累积下渗量 [cm]
  ponding_time::T = NaN  # 积水时间
end

function infiltration_GreenAmpt!(ga::GreenAmpt, P::T, dt::T) where T
  """
  P: 降雨强度 [cm/h]
  dt: 时间步长 [s]
  """
  Δθ = ga.θs - ga.θi
  dt_h = dt / 3600.0

  if isnan(ga.ponding_time)
    # 积水前: f = P
    if P <= ga.Ksat * (1 + ga.ψf * Δθ / ga.F)
      f = P
      ga.F += f * dt_h
    else
      # 发生积水
      ga.ponding_time = t
      ga.F = ga.Ksat * ga.ψf * Δθ / (P - ga.Ksat)
    end
  else
    # 积水后: Green-Ampt 方程
    # f = Ksat * (1 + ψf * Δθ / F)
    f = ga.Ksat * (1 + ga.ψf * Δθ / ga.F)
    ga.F += f * dt_h
  end

  return f  # 实际下渗率
end
```

## 方案 3: 完整的地表-土壤耦合模块

整合到 `Soil` 结构中：

### 新增字段 (在 `Soil.jl`):
```julia
# 下渗模块
P::FT = 0.0              # 降雨强度 [cm/h]
infiltration::FT = 0.0   # 实际下渗率 [cm/h]
ponding_depth::FT = 0.0  # 积水深度 [cm]
runoff::FT = 0.0         # 地表径流 [cm/h]
```

### 新函数 (`Infiltration.jl`):
```julia
function update_infiltration!(soil::Soil, P::Float64, dt::Float64)
  """
  完整的降雨-下渗-径流过程
  """

  # 1. 计算下渗能力
  infiltration_capacity = calc_infiltration_capacity(soil)

  # 2. 水量平衡
  if soil.ponding_depth > 0
    # 有积水：全速下渗
    actual_inf = min(infiltration_capacity,
                     (P + soil.ponding_depth/(dt/3600)))
    soil.ponding_depth += (P - actual_inf) * dt/3600
  else
    # 无积水
    actual_inf = min(P, infiltration_capacity)

    # 检查是否会产生积水
    if P > infiltration_capacity
      excess = (P - infiltration_capacity) * dt/3600
      soil.ponding_depth = excess
    end
  end

  # 3. 超过最大积水深度则产生径流
  ponding_max = 0.5  # cm
  if soil.ponding_depth > ponding_max
    soil.runoff = (soil.ponding_depth - ponding_max) / (dt/3600)
    soil.ponding_depth = ponding_max
  end

  # 4. 更新边界条件
  soil.Q0 = -actual_inf
  soil.infiltration = actual_inf

  return actual_inf
end
```

## 方案 4: 冻土下渗模块

如果考虑季节冻土：

```julia
function infiltration_frozen!(soil::Soil, P::Float64)
  """
  考虑冻结层的下渗
  """

  # 检查冻结层
  frozen_layers = soil.Tsoil .< 0.0

  if any(frozen_layers[1:3])  # 表层冻结
    # 冻土下渗能力折减
    ice_fraction = calc_ice_content(soil)
    reduction_factor = (1 - ice_fraction)^3  # 经验公式

    K_frozen = soil.K[1] * reduction_factor
    infiltration = min(P, K_frozen)

    # 超渗产流
    soil.runoff = P - infiltration
    soil.Q0 = -infiltration
  else
    # 正常下渗
    infiltration = P
    soil.Q0 = -infiltration
  end
end
```

## 推荐实施方案

### 分阶段实施:

#### 第一阶段（最简单，立即可用）
- [ ] 在 `Equation_Richards.jl` 中添加 `infiltration_boundary!` 函数
- [ ] 修改 `RichardsEquation` 调用它来自适应切换边界
- [ ] 在 `Soil` 中添加 `P`, `infiltration`, `runoff` 字段
- [ ] 编写基础测试用例

#### 第二阶段（增强）
- [ ] 创建新文件 `src/SoilMoisture/Infiltration.jl`
- [ ] 实现 Green-Ampt 作为快速替代
- [ ] 添加积水深度追踪
- [ ] 完善地表-土壤耦合
- [ ] 添加更多测试用例

#### 第三阶段（可选）
- [ ] 冻土下渗
- [ ] 优先流（macropore flow）
- [ ] 地表糙度影响
- [ ] 性能优化和基准测试

### 测试用例建议

创建 `test/test-infiltration.jl`:

```julia
using SoilDifferentialEquations, Test

@testset "Infiltration Module" begin

  # Test 1: 恒定降雨，P < Ksat
  @testset "Constant rainfall below capacity" begin
    soil = init_test_soil()
    P = 2.0  # cm/h，中等强度
    dt = 60.0  # seconds

    method = infiltration_boundary!(soil, P, dt)
    @test method == "Q0"
    @test soil.Q0 ≈ -P
    @test soil.runoff == 0.0
  end

  # Test 2: 强降雨超渗，P > Ksat
  @testset "Heavy rainfall exceeding capacity" begin
    soil = init_test_soil()
    P = 10.0  # cm/h，超过 Ksat
    dt = 60.0

    # 多个时间步
    for i = 1:10
      method = infiltration_boundary!(soil, P, dt)
    end

    @test soil.ponding_depth > 0  # 应该有积水
    @test soil.infiltration < P   # 实际下渗 < 降雨
  end

  # Test 3: 间歇降雨
  @testset "Intermittent rainfall" begin
    soil = init_test_soil()
    dt = 60.0

    # 1小时强降雨
    for i = 1:60
      infiltration_boundary!(soil, 5.0, dt)
    end
    ponding_after_rain = soil.ponding_depth

    # 2小时无降雨
    for i = 1:120
      infiltration_boundary!(soil, 0.0, dt)
    end

    @test soil.ponding_depth < ponding_after_rain  # 积水应减少
  end

  # Test 4: 质量守恒检验
  @testset "Mass balance check" begin
    soil = init_test_soil()
    P = 3.0
    dt = 60.0
    ntim = 100

    total_rain = 0.0
    total_infiltration = 0.0
    total_runoff = 0.0

    for i = 1:ntim
      infiltration_boundary!(soil, P, dt)
      total_rain += P * dt/3600
      total_infiltration += soil.infiltration * dt/3600
      total_runoff += soil.runoff * dt/3600
    end

    final_ponding = soil.ponding_depth
    mass_balance = total_rain - total_infiltration - total_runoff - final_ponding

    @test abs(mass_balance) < 0.01  # 质量守恒误差 < 0.01 cm
  end

end
```

### 关键设计考虑

1. **边界条件切换**
   - 需要平滑过渡，避免数值震荡
   - 考虑滞后效应（积水消退后不立即切换）

2. **时间步长限制**
   - 下渗过程可能需要更小的时间步长
   - 考虑自适应时间步长

3. **物理约束**
   - θ不能超过θ_sat
   - 积水深度不能为负
   - 所有通量守恒

4. **性能优化**
   - Green-Ampt 比完整 Richards 快很多
   - 可以根据情况选择不同复杂度的模型

## 参考文献

1. Infiltration Theory:
   - Green, W.H. and Ampt, G.A. (1911)
   - Philip, J.R. (1957) - 两项入渗方程

2. Numerical Methods:
   - Celia et al. (1990) - Mixed form Richards equation
   - Kavetski et al. (2001) - Adaptive time stepping

3. Frozen Soil:
   - Zhao et al. (1997) - 冻土水热耦合
   - Hansson et al. (2004) - 相变处理