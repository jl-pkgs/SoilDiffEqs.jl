# SoilDiffEqs.jl Agent Guidelines

SoilDiffEqs.jl 是一个用于求解土壤水热运动微分方程的 Julia 软件包，实现了 Richards 方程、土壤热传导方程以及地下水动态模拟。

> **代码编写准则**：

**(a) 代码编写遵循Linux极简主义哲学，同时符合代码排版规范**

Linux极简主义哲学的核心原则：

1. 一个程序只做一件事，并把它做好
2. 简洁至上
3. 不要重复造轮子

应用到Julia代码中：

1. 函数应该短小精悍，只做一件事
2. 避免不必要的抽象和过度设计
3. 清晰的命名，自解释的代码
4. 减少重复代码
5. 使用标准库和现有工具

**(b) 代码责任制**

谁修改谁负责，修改后的代码，一定要测试。

改完代码不测试，是不负责任的流氓。


## julia独有的语法

```julia
# 不推荐
zs_obs_orgin = config.zs_obs_orgin
zs_obs = config.zs_obs

# 推荐
(; zs_obs_orgin, zs_obs) = config
```

---

## 1. 项目概述

### 1.1 核心功能
- **土壤水分模拟**: 基于 Richards 方程求解土壤水分的一维垂直运动
- **土壤温度模拟**: 基于热传导方程计算土壤温度剖面
- **地下水模拟**: 模拟地下水位动态及其与土壤水分的相互作用
- **植被-土壤耦合**: 考虑根系吸水和土壤水分限制因子

### 1.2 技术栈
- **语言**: Julia (v1.10+)
- **核心依赖**:
  - `DifferentialEquations.jl` / `OrdinaryDiffEqTsit5`: ODE 求解器
  - `Parameters.jl`: 参数结构化
  - `OffsetArrays.jl`: 偏移数组索引
  - `StructArrays.jl`: 结构化数组
  - `ModelParams.jl`: 参数优化 (SCE-UA 算法)

### 1.3 物理基础
- **Richards 方程**: 描述土壤水分非饱和流动的控制方程
- **达西定律**: 土壤水流通量计算
- **土壤水分特征曲线**: 支持 van Genuchten 和 Campbell 两种参数化方案
- **热传导方程**: 土壤温度传输

---

## 2. 项目结构

```
SoilDiffEqs.jl/
├── Project.toml              # Julia 项目配置
├── Manifest.toml             # 依赖锁定文件
├── Artifacts.toml            # 数据 artifacts 配置
├── codecov.yml               # 代码覆盖率配置
├── src/
│   ├── SoilDifferentialEquations.jl    # 主模块入口
│   ├── GlobalOptions.jl                # 全局选项配置
│   ├── Soil.jl                         # 核心 Soil 数据结构
│   ├── SoilParam.jl                    # 土壤参数定义
│   ├── Config.jl                       # 配置管理，支持 YAML 配置文件（新增）
│   ├── soil_texture.jl                 # USDA 土壤质地分类
│   ├── tridiagonal_solver.jl           # 三对角矩阵求解器
│   ├── solver.jl                       # 自定义 ODE solver (实验性)
│   ├── case_Bonan2019.jl               # Bonan 2019 测试案例
│   ├── SPAC.jl                         # 陆面过程综合模拟器
│   │
│   ├── SoilMoisture/                   # 土壤水分模块
│   │   ├── SoilMoisture.jl             # 子模块入口
│   │   ├── Equation_Richards.jl        # Richards 方程核心实现
│   │   ├── Equation_Zeng2009.jl        # Zeng 2009 方案
│   │   ├── Retention.jl                # 土壤水分特征曲线接口
│   │   ├── Retention_van_Genuchten.jl  # van Genuchten 模型
│   │   ├── Retention_Campbell.jl       # Campbell 模型
│   │   ├── Equilibrium.jl              # 平衡土壤含水量计算
│   │   ├── soil_moisture_Bonan.jl      # Bonan 隐式求解方案
│   │   ├── soil_moisture_Bonan_Q0.jl   # Bonan 方案 (通量边界)
│   │   ├── soil_moisture_Zeng2009.jl   # Zeng 2009 求解方案
│   │   ├── soil_moisture_BEPS.jl       # BEPS 模型方案
│   │   ├── Solve_SM.jl                 # 土壤 moisture 求解器接口
│   │   └── clm5_utils.jl               # CLM5 相关工具函数
│   │
│   ├── SoilTemperature/                # 土壤温度模块
│   │   ├── soil_properties_thermal.jl  # 热力参数计算
│   │   ├── soil_temperature.jl         # 温度求解 (Dirichlet 边界)
│   │   ├── soil_temperature_F0.jl      # 温度求解 (Neumann 边界)
│   │   ├── EquationTsoil.jl            # 热传导方程实现
│   │   └── Solve_Tsoil.jl              # 温度求解器接口
│   │
│   ├── GroundWater/                    # 地下水模块
│   │   ├── GroundWater.jl              # 子模块入口
│   │   ├── GW_Update_ZWT.jl            # 地下水位更新
│   │   ├── GW_Correctθ.jl              # 含水量修正
│   │   └── specific_yield.jl           # 给水度计算
│   │
│   └── Vegetation/                     # 植被模块
│       ├── root_fraction.jl            # 根系分布计算
│       └── soil_moisture_constraint.jl # 土壤水分限制因子
│
├── ext/
│   └── SoilDiffEqExt.jl                # 扩展模块 (ODE 求解器集成)
│
├── test/                               # 测试目录
│   ├── runtests.jl                     # 测试入口
│   ├── test-soil_moisture.jl           # 土壤 moisture 测试
│   ├── test-soil_moisture_Q0.jl        # 通量边界测试
│   ├── test-soil_temperature.jl        # 温度求解测试
│   ├── test-solve_SM.jl                # 求解器接口测试
│   ├── GW/                             # 地下水测试
│   └── SM_uscrn/                       # USCRN 数据测试
│
├── examples/                           # 示例代码
│   ├── SM_uscrn/                       # USCRN 站点模拟
│   ├── SM_China/                       # 中国站点模拟（已重构，支持 YAML 配置）
│   │   ├── run_config.jl               # 基于配置的主运行脚本
│   │   ├── config.yaml                 # YAML 配置文件
│   │   ├── src/main_plot.jl            # 绘图功能
│   │   └── src/main_config.jl          # 配置辅助函数
│   ├── Tsoil_ex01_CUG/                 # 温度模拟示例
│   └── common/SoilConfig.jl            # 配置工具
│
└── docs/                               # 文档
    ├── 模型手册.md                      # 详细使用手册
    ├── SM_solver.md                    # 求解器文档
    └── ...
```

---

## 3. 构建与测试

### 3.1 运行测试

```bash
# 运行全部测试
julia --project -e "using Pkg; Pkg.test()"

# 运行单个测试文件
julia --project test/test-soil_moisture.jl

# 交互式运行测试
julia --project
] test
```

### 3.2 开发环境设置

```julia
# 激活项目环境
using Pkg
Pkg.activate(".")

# 安装依赖
Pkg.instantiate()

# 加载模块
using SoilDifferentialEquations
```

### 3.3 CI/CD

- **平台**: GitHub Actions (`.github/workflows/CI.yml`)
- **Julia 版本**: v1 (最新稳定版)
- **操作系统**: Ubuntu Linux
- **代码覆盖率**: Codecov 集成

---

## 4. 代码风格指南

### 4.1 命名规范

- **函数**: 使用 `snake_case`，突变函数以 `!` 结尾
  - 示例: `cal_Q!`, `soil_Updateθ!`, `soil_moisture!`
- **类型/结构体**: 使用 `CamelCase`
  - 示例: `Soil`, `SoilParam`, `VanGenuchten`, `Campbell`
- **变量**: 混合使用 `snake_case` 和数学符号
  - 希腊字母: `θ` (含水量), `ψ` (水势), `Δz` (层厚), `κ` (热导率)
  - 下标符号: `K₊ₕ` (界面导水率), `Δz₊ₕ` (界面间距), `z₊ₕ` (界面深度)
  - 常规变量: `ibeg` (起始活跃层索引), `zwt` (地下水位深度)

### 4.2 单位约定

| 变量       | 符号 | 单位     | 说明             |
| ---------- | ---- | -------- | ---------------- |
| 体积含水量 | θ    | [m³ m⁻³] | 土壤水分体积比   |
| 土壤水势   | ψ    | [cm]     | 负值表示吸力     |
| 水流通量   | Q    | [cm h⁻¹] | 向下为负         |
| 深度       | z    | [m]      | 向下为负         |
| 深度       | z_cm | [cm]     | 厘米单位版本     |
| 时间步长   | dt   | [s]      | 计算时转换为小时 |
| 热通量     | F    | [W m⁻²]  | 向下为负         |
| 温度       | T    | [°C]     | 摄氏度           |

### 4.3 类型系统

- 使用参数化类型保证性能和精度控制
```julia
function cal_Q!(soil::Soil{T}; ...) where {T<:Real}
```
- 使用多重派发处理不同土壤参数模型
```julia
function clamp_θ!(soil::Soil{T,VanGenuchten{T}}, ...) where {T<:Real}
function clamp_θ!(soil::Soil{T,Campbell{T}}, ...) where {T<:Real}
```

### 4.4 性能优化

- 使用 `@inbounds` 优化热循环
- 避免在内部循环中分配内存 (预分配数组在 `Soil` 结构体中)
- 使用 `(; var1, var2) = struct` 语法解包结构体
- 使用 `OffsetArrays` 实现方便的 0-based 索引

### 4.5 注释与文档

- 文档字符串使用 Markdown 格式
- 注释可使用中文或英文
- 复杂算法需要说明物理意义和数学公式

---

## 5. 核心数据结构

### 5.1 Soil 结构体

`Soil` 是核心数据结构，包含以下主要字段:

```julia
@with_kw_noshow mutable struct Soil{FT,P<:AbstractSoilParam{FT}}
  # 网格配置
  N::Int                        # 土壤层数
  ibeg::Int                     # 第一个活跃层索引
  z::OffsetVector{FT}           # 层中心深度 [m], 向下为负
  Δz::Vector{FT}                # 层厚度 [m]
  z₊ₕ::Vector{FT}               # 界面深度 [m]
  
  # 水分变量
  θ::Vector{FT}                 # 体积含水量 [m³/m³]
  ψ::Vector{FT}                 # 土壤水势 [cm]
  K::Vector{FT}                 # 水力传导度 [cm/h]
  K₊ₕ::Vector{FT}               # 界面传导度 [cm/h]
  Q::Vector{FT}                 # 水流通量 [cm/h]
  sink::Vector{FT}              # 源汇项 (蒸散发) [cm/time]
  
  # 边界条件
  θ0::FT                        # 地表含水量
  ψ0::FT                        # 地表水势
  Q0::FT                        # 地表通量
  
  # 地下水
  zwt::FT                       # 地下水位深度 [m]
  wa::FT                        # 含水层储水量 [mm]
  
  # 温度变量
  Tsoil::Vector{FT}             # 土壤温度 [°C]
  Tsurf::FT                     # 地表温度 [°C]
  F0::FT                        # 地表热通量 [W/m²]
  
  # 植被
  f_root::Vector{FT}            # 根系分布比例
  β::FT                         # 土壤水分限制因子
  
  # 参数
  param::SoilParam{FT,P}        # 土壤水力+热力参数
end
```

### 5.2 Config 结构体（新增）

用于 YAML 配置文件管理:

```julia
@with_kw struct Config
  # data: 数据相关配置
  file::String              # 数据文件路径
  col_time::Int = 1         # 时间列索引
  col_obs_start::Int = 2    # 观测数据起始列
  scale_factor::Float64 = 1.0
  zs_obs_orgin::Vector{Float64}  # 原始观测深度 [cm]
  zs_obs::Vector{Float64}        # 插值后观测深度 [cm]
  z_bound_top::Float64 = 10.0    # 顶层边界层深度 [cm]

  # model: 模型参数配置
  soil_type::Int = 7
  same_layer::Bool = false
  method_retention::String = "van_Genuchten"  # 持水曲线方法
  method_solve::String = "Bonan"              # 求解方法: "Bonan" 或 "ODE"
  dt::Float64 = 3600.0
  zs_center::Vector{Float64}  # 模拟层中心深度 [cm]

  # optimization: 优化配置
  optim::Bool = false
  maxn::Int = 100
  objective::String = "NSE"   # 目标函数: "NSE" 或 "KGE"

  # output: 输出配置
  plot_file::String = ""
end
```

使用 `load_config(cfg_file)` 从 YAML 加载配置。

### 5.3 网格系统 (交错网格)

- **状态变量** (θ, ψ, T): 定义在层中心 `z[i]`
- **通量变量** (Q, K, F): 定义在层界面 `z₊ₕ[i]`
- 界面位于层 i 和层 i+1 之间

```
      z[1]        z[2]              z[N]
        |           |                  |
   =====|===========|==================|===== 地表 (z=0)
        |     ↑     |                  |
       θ[1]  K₊ₕ[1] θ[2]  ...         θ[N]
        |     |     |
        | --- | --- ||  |z₊ₕ[N] (底部)
```

### 5.4 土壤参数模型

支持两种土壤水分特征曲线模型:

1. **van Genuchten 模型** (推荐)
   - 参数: θ_sat, θ_res, Ksat, α, n, m
   - 公式: `θ(ψ) = θ_r + (θ_s - θ_r) / [1 + |αψ|ⁿ]ᵐ`

2. **Campbell 模型** (简化)
   - 参数: θ_sat, ψ_sat, Ksat, b
   - 公式: `θ(ψ) = θ_s * (ψ/ψ_sat)^(-1/b)`

---

## 6. 主要求解器

### 6.1 土壤水分求解

| 函数                      | 边界条件  | 算法       | 用途                          |
| ------------------------- | --------- | ---------- | ----------------------------- |
| `soil_moisture!`          | ψ0 (水势) | Bonan 隐式 | 积水/灌溉条件                 |
| `soil_moisture_Q0!`       | Q0 (通量) | Bonan 隐式 | 降雨入渗条件                  |
| `soil_moisture_Zeng2009!` | 两者      | Zeng 2009  | 通用求解                      |
| `RichardsEquation`        | 两者      | ODE        | 与 DifferentialEquations 集成 |

### 6.2 土壤温度求解

| 函数                   | 边界条件          | 说明           |
| ---------------------- | ----------------- | -------------- |
| `soil_temperature!`    | Tsurf (Dirichlet) | 已知地表温度   |
| `soil_temperature_F0!` | F0 (Neumann)      | 已知地表热通量 |

---

## 7. 常见任务指南

### 7.1 创建土壤对象

```julia
# 从层厚度创建
Δz = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]  # 6层土壤 [m]
soil = Soil(Δz; method_retention="van_Genuchten")

# 直接构造
soil = Soil{Float64}(; N=10, method_retention="Campbell")
```

### 7.2 设置土壤参数

```julia
# van Genuchten 参数
soil.param.θ_sat .= 0.40    # 饱和含水量
soil.param.θ_res .= 0.078   # 残余含水量
soil.param.Ksat .= 1.04     # 饱和导水率 [cm/h]
soil.param.α .= 0.036       # 进气值倒数 [cm⁻¹]
soil.param.n .= 1.56        # 孔隙分布指数
```

### 7.3 运行单步模拟

```julia
# 水势边界 (ψ0 = 0 表示地表积水)
soil.ψ0 = 0.0
soil.sink .= 0.01  # 蒸发速率 [cm/h]
soil_moisture!(soil, soil.sink, soil.ψ0)

# 通量边界 (降雨入渗)
Q0 = -2.0  # 入渗速率 2 cm/h (负值表示向下)
soil_moisture_Q0!(soil, soil.sink, Q0)
```

### 7.4 使用 ODE 求解器

```julia
using OrdinaryDiffEqTsit5

u0 = soil.θ[soil.ibeg:soil.N]  # 初始状态
tspan = (0.0, 3600.0)          # 时间范围 [s]
prob = ODEProblem(RichardsEquation, u0, tspan, soil)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
```

### 7.5 使用 YAML 配置运行（新增）

```julia
# 加载配置
config = load_config("examples/SM_China/config.yaml")

# 运行完整模拟（包含数据加载、模拟、优化、绘图）
include("examples/SM_China/run_config.jl")
```

配置文件示例 (`config.yaml`):
```yaml
data:
  file: "theta/theta_647.csv"
  col_time: 1
  col_obs_start: 2
  scale_factor: 0.01
  zs_obs_orgin: [5, 10, 20, 50, 100]
  zs_obs: [10, 20, 50, 100]      # 顶层10cm作为边界条件
  z_bound_top: 10.0

model:
  soil_type: 7
  same_layer: false
  method_retention: "van_Genuchten"
  method_solve: "Bonan"          # 或 "ODE"
  dt: 3600.0
  zs_center: [0, -5, -15, -30, -60, -100, -150, -200]

optimization:
  enable: true
  maxn: 100
  objective: "NSE"               # 或 "KGE"

output:
  plot_file: "images/plot.png"
```

---

## 8. 注意事项

### 8.1 不应修改的内容

- **Unicode 数学符号**: 保留 `θ`, `ψ`, `Δz`, `K₊ₕ` 等，不要改为 ASCII 名称
- **`@inbounds` 注解**: 不要未经验证就移除
- **`Soil` 结构体布局**: 修改前需考虑内存布局影响

### 8.2 数值稳定性

- 时间步长 `dt` 建议 300-3600 秒 (5分钟到1小时)
- 降雨入渗模拟需要更小步长 (60-600 秒)
- 含水量会被自动限制在物理有效范围内

### 8.3 质量守恒检查

```julia
# 检查水量平衡误差
info = error_SM(soil)
# info.bias: 绝对误差 [cm]
# info.perc: 相对误差 [%]
```

---

## 9. 扩展模块 (Weak Dependencies)

项目使用 Julia 的扩展系统 (package extensions):

- **SoilDiffEqExt**: 当 `OrdinaryDiffEqTsit5` 加载时自动激活
  - 提供 `_ODEProblem` 和 `_solve` 接口
  - 桥接自定义类型与 DifferentialEquations.jl

---

## 10. 参考文档

- 详细模型手册: `docs/模型手册.md`
- 求解器文档: `docs/SM_solver.md`
- 土壤特征曲线: `docs/Retention.md`
- 地下水模块: `docs/ch09_groundwater.md`

---

**最后更新**: 2026-02-02
