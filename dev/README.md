# 地表积水处理方案

本文件夹包含三种地表积水处理方案的实现代码，用于模拟地表积水的产流过程。

## 📁 文件说明

| 文件 | 说明 |
|------|------|
| `ponding_schemeA.jl` | 方案A：基于蓄水容量的方案 |
| `ponding_schemeB.jl` | 方案B：基于坡面汇流的方案（曼宁公式） |
| `ponding_schemeC.jl` | 方案C：混合方案（推荐） |
| `test_ponding_schemes.jl` | 测试和示例代码 |
| `README.md` | 本说明文档 |

## 🎯 方案对比

### 方案A：基于蓄水容量

**原理**：定义最大地表蓄水深度，超过部分立即产流（蓄满产流）

**适用场景**：
- ✅ 平坦地区
- ✅ 农田
- ✅ 城市下垫面

**优点**：
- 物理意义明确
- 计算简单高效
- 参数少（仅1个：`h_pond_max`）

**参数**：
```julia
h_pond_max = 0.5  # 最大蓄水深度 [cm]
```

### 方案B：基于坡面汇流

**原理**：使用曼宁公式计算径流速度，考虑坡度、糙度、积水深度的非线性影响

**适用场景**：
- ✅ 坡地
- ✅ 山区
- ✅ 自然流域

**优点**：
- 物理过程清晰
- 考虑地形影响
- 积水深度影响强（h^5/3）

**参数**：
```julia
S0 = 0.05          # 坡度 [-]
n_manning = 0.15   # 曼宁粗糙系数 [s/m^(1/3)]
L_char = 50.0      # 特征汇流长度 [m]
h_pond_min = 0.1   # 起始产流深度 [cm]
```

### 方案C：混合方案 ⭐ 推荐

**原理**：三阶段机制
1. h < h_pond_min：微地形蓄水，无径流
2. h_pond_min < h < h_pond_max：坡面汇流（曼宁公式）
3. h > h_pond_max：超蓄产流

**适用场景**：
- ✅ 通用场景
- ✅ 需要综合考虑多种过程
- ✅ 有观测数据可率定

**优点**：
- 物理机制完整
- 适应性强
- 参数可率定

**参数**：
```julia
h_pond_min = 0.1   # 微地形蓄水容量 [cm]
h_pond_max = 2.0   # 最大允许积水深度 [cm]
S0 = 0.01          # 坡度 [-]
n_manning = 0.15   # 曼宁系数 [s/m^(1/3)]
L_char = 100.0     # 汇流特征长度 [m]
```

## 🚀 快速开始

### 1. 加载代码

```julia
# 在Julia REPL中
cd("path/to/SoilDiffEqs.jl")

# 加载单个方案
include("dev/ponding_schemeA.jl")

# 或加载所有方案
include("dev/ponding_schemeA.jl")
include("dev/ponding_schemeB.jl")
include("dev/ponding_schemeC.jl")
```

### 2. 基本使用

```julia
# 方案A：简单场景
h_current = 0.3   # 当前积水深度 [cm]
P = 10.0          # 降水强度 [cm/h]
inf = 3.0         # 下渗速率 [cm/h]
dt = 0.1          # 时间步长 [h]

runoff, h_new = ponding_schemeA(h_current, P, inf, dt; h_pond_max=0.5)
println("径流: $runoff cm, 新积水: $h_new cm")

# 方案B：坡地场景
runoff, h_new = ponding_schemeB(h_current, P, inf, dt;
                                S0=0.05, n_manning=0.15, L_char=50.0)

# 方案C：通用场景（推荐）
runoff, h_new = ponding_schemeC(h_current, P, inf, dt;
                                h_pond_min=0.1, h_pond_max=2.0,
                                S0=0.01, n_manning=0.15, L_char=100.0)
```

### 3. 运行测试

```julia
include("dev/test_ponding_schemes.jl")

# 运行所有测试
test_all_schemes()

# 或运行单个测试
test_schemeA()
test_schemeB()
test_schemeC()
compare_schemes()
simulate_rainfall_event()
```

## 📊 使用示例

### 示例1：平坦农田 - 间歇性降雨

```julia
include("dev/ponding_schemeA.jl")

# 初始化
h_pond = 0.0
h_pond_max = 0.5  # [cm]
dt = 1/6  # 10分钟 [h]

# 模拟3小时降雨
P_series = [10.0, 15.0, 20.0, 5.0, 2.0, 0.0] .* ones(3)  # 18步
inf_series = [5.0, 4.0, 3.0, 2.5, 2.0, 1.5] .* ones(3)

for (P, inf) in zip(P_series, inf_series)
    runoff, h_pond = ponding_schemeA(h_pond, P, inf, dt; h_pond_max=h_pond_max)
    println("P=$P, inf=$inf → h=$h_pond cm, runoff=$runoff cm")
end
```

### 示例2：坡地 - 单次暴雨

```julia
include("dev/ponding_schemeB.jl")

# 坡地参数
S0 = 0.08          # 8% 坡度
n_manning = 0.18   # 农田糙度
L_char = 40.0      # 40m坡长

h_pond = 0.0
dt = 0.05  # 3分钟 [h]

# 暴雨30分钟
for i in 1:10
    P = 30.0  # [cm/h]
    inf = max(5.0 - i*0.3, 1.5)  # 下渗递减

    runoff, h_pond = ponding_schemeB(h_pond, P, inf, dt;
                                     S0=S0, n_manning=n_manning, L_char=L_char)
    println("t=$(i*3)min: h=$(round(h_pond,digits=2)) cm, runoff=$(round(runoff,digits=3)) cm")
end
```

### 示例3：与Soil对象配合使用

```julia
using Parameters  # 如果有的话

include("dev/ponding_schemeC.jl")
# include("src/Soil.jl")  # 加载Soil结构

# 假设已有soil对象
# soil = Soil(...)

# 在时间循环中调用
h_pond = 0.0

for t in 1:n_steps
    P = rainfall[t]      # 降水强度 [cm/h]
    inf = infiltration[t]  # 下渗速率 [cm/h]

    # 使用!版本（自动从soil.dt获取时间步长）
    runoff, h_pond = ponding_schemeC!(soil, h_pond, P, inf;
                                      h_pond_min=0.1, h_pond_max=2.0,
                                      S0=0.02, n_manning=0.15, L_char=80.0)

    # 更新其他变量...
end
```

## 🔧 参数选择指南

### 方案A参数

| 下垫面类型 | h_pond_max [cm] |
|------------|-----------------|
| 城市不透水面 | 0.1 - 0.3 |
| 平坦农田 | 0.3 - 0.5 |
| 有起伏的地形 | 0.5 - 1.0 |

### 方案B参数

#### 曼宁粗糙系数 n [s/m^(1/3)]

| 地表类型 | n值范围 |
|----------|---------|
| 裸地/光滑表面 | 0.01 - 0.05 |
| 农田 | 0.10 - 0.20 |
| 草地 | 0.15 - 0.40 |
| 森林/灌木 | 0.30 - 0.80 |

#### 坡度 S₀ [-]

| 地形 | S₀ | 百分比 |
|------|-----|--------|
| 平坦 | 0.001 - 0.01 | 0.1% - 1% |
| 缓坡 | 0.01 - 0.05 | 1% - 5% |
| 中坡 | 0.05 - 0.15 | 5% - 15% |
| 陡坡 | 0.15 - 0.30 | 15% - 30% |

#### 特征长度 L [m]（单站点模拟）

| 场景 | L值 |
|------|-----|
| 农田/草地 | 30 - 50 |
| 森林 | 50 - 100 |
| 城市 | 10 - 30 |
| 快速测试 | 50 |

### 方案C参数

综合方案A和B的参数选择指南。

## 📖 理论文档

详细的理论推导和公式说明见：
```
docs/下渗/Infiltration_final.typ
```

包含：
- 曼宁公式推导
- q/L的物理意义
- 典型排水速率参考值
- 单站点L参数取值指南

## 🔬 与Richards方程的耦合

所有方案计算的积水深度 `h_pond` 可以直接用于Richards方程的边界条件：

```julia
# 最大下渗能力（1.1节）
f_potential = -K₁ * (2*(h_pond - ψ₁)/Δz₁ + 1)
```

**耦合机制**：
- h_pond ↑ → 水头梯度 ↑ → 下渗能力 ↑
- 下渗量 ↑ → 积水消耗 ↑ → h_pond ↓
- 形成负反馈调节

## ⚠️ 注意事项

1. **单位一致性**
   - 输入输出：[cm] 和 [h]
   - 内部计算：部分使用 [m] 和 [s]
   - 已自动处理转换

2. **时间步长**
   - 推荐：dt ≤ 0.1 h（6分钟）
   - 暴雨事件：dt ≤ 0.05 h（3分钟）

3. **数值稳定性**
   - 方案B/C在h很小时已做保护处理
   - 避免负积水深度（已自动处理）

4. **参数率定**
   - 优先使用观测数据率定
   - 无观测时使用典型值
   - 做敏感性分析

## 🔜 后续工作

测试通过后，可将代码集成到主模型：

1. 在 `Soil` 结构中添加积水相关字段：
```julia
h_pond::FT = FT(0.0)  # 地表积水深度 [cm]
```

2. 创建统一接口函数

3. 添加到时间积分循环中

4. 编写单元测试

## 📝 更新日志

- 2025-11-05: 初始版本，实现三种方案

## 📧 联系方式

如有问题，请参考：
- 理论文档：`docs/下渗/Infiltration_final.typ`
- 代码注释：各方案文件的文档字符串

---

**祝使用愉快！** 🎉
