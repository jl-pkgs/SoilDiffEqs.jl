# SoilDiffEqs.jl Agent Guidelines

> **最后更新**: 2026-02-05, commit: b53d0e0 

SoilDiffEqs.jl 是一个用于求解土壤水热运动微分方程的 Julia 软件包，实现了 Richards 方程、土壤热传导方程以及地下水动态模拟。

> **AGENTS.md**更新方法
读取最近3次的commit，翻阅相关代码，更新AGENTS.md

代码编写准则见：`.github/Coder.md`。这是代码编写的第一宪法，请务必遵循。


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
│   ├── Config.jl                       # 配置管理，支持 YAML 配置文件
│   ├── Model_Soil.jl                   # 基于配置的运行接口（setup, Soil_predict, Soil_main 等）
│   ├── soil_texture.jl                 # USDA 土壤质地分类（枚举类型、中英文支持）
│   ├── tridiagonal_solver.jl           # 三对角矩阵求解器
│   ├── solver.jl                       # 自定义 ODE solver (实验性)
│   ├── ultilize.jl                     # 工具函数（深度插值、索引查找）
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
│   ├── SM/                             # 土壤水分模拟示例（YAML 配置驱动）
│   │   ├── case_SM_Bonan_China.jl      # Bonan求解器中国站点示例
│   │   ├── case_SM_Bonan_China.yaml    # Bonan求解器中国站点配置
│   │   ├── case_SM_Bonan_uscrn.jl      # Bonan求解器USCRN示例
│   │   ├── case_SM_Bonan_uscrn.yaml    # Bonan求解器USCRN配置
│   │   ├── BEPS/                       # BEPS求解器示例
│   │   │   ├── case_SM_BEPS_uscrn.jl
│   │   │   └── case_SM_BEPS_uscrn.yaml
│   │   ├── data/                       # 数据目录
│   │   ├── images/                     # 输出图表目录
│   │   └── output/                     # 优化结果输出目录
│   ├── Tsoil/                          # 土壤温度模拟示例
│   │   ├── case_Tsoil_CUG.jl           # CUG 站点温度模拟
│   │   ├── case_Tsoil_CUG.yaml         # CUG 配置文件
│   │   ├── case_Tsoil_isusm.jl         # ISUSM 温度模拟
│   │   ├── case_Tsoil_isusm.yaml       # ISUSM 配置文件
│   │   ├── case_Tsoil_uscrn_03.jl      # USCRN 温度模拟
│   │   └── case_Tsoil_uscrn_03.yaml    # USCRN 温度配置文件
│   ├── main_plot.jl                    # 通用绘图功能
│   └── common/SoilConfig.jl            # 配置工具（旧版）
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
@with_kw mutable struct Config
  ## 模型类型: "SM" | "Tsoil"
  model_type::String = "SM"

  ## data: 数据相关配置
  file::String = ""              # 数据文件路径
  fileConfig::String = ""        # 配置文件路径
  col_time::Int = 1              # 时间列索引
  col_obs_start::Int = 2         # 观测数据起始列
  scale_factor::Float64 = 1.0
  zs_obs_orgin::Vector{Float64}  # 原始观测深度 [cm]，Tsoil可选
  zs_obs::Vector{Float64}        # 插值后观测深度 [cm]，Tsoil可选
  z_bound_top::Float64 = 10.0    # 顶层边界层深度 [cm]
  itop::Int = 1                  # [derived] 上边界层在观测数据中的索引

  ## model: 模型参数配置
  soil_type::Int = 7
  same_layer::Bool = false
  method_retention::String = "van_Genuchten"  # 持水曲线方法（SM特有）
  method_solve::String = "Bonan"              # 求解方法: "Bonan" | "ODE" | "BEPS"
  dt::Float64 = 3600.0
  zs_center::Vector{Float64}     # 模拟层中心深度 [cm]
  grid::NamedTuple               # [derived] 预计算的网格信息 (z, z₊ₕ, Δz, Δz₊ₕ, N, ibeg, itop)
  
  # 自定义土壤参数（可选，用于覆盖标准参数）- SM 特有
  soil_params::Union{Dict{String,Float64},Nothing} = nothing

  ## optimization: 优化配置
  optim::Bool = false
  maxn::Int = 100
  objective::String = "NSE"   # 目标函数: "NSE" 或 "KGE"
  of_fun::Function = of_NSE   # 目标函数实现（根据 objective 自动设置）

  # output: 输出配置
  plot_file::String = ""

  # internal
  config_file::String = ""  # 配置文件路径（用于日志命名）
end
```

使用 `load_config(fileConfig)` 从 YAML 加载配置。配置文件中 `optimization.enable` 映射到 `optim` 字段，`optimization.objective` 决定 `of_fun` 的值（NSE→`of_NSE`，KGE→`of_KGE`）。

**Config 网格计算**: `load_config` 会自动计算并填充 `grid` 字段，包含：
- `z`, `z₊ₕ`, `Δz`, `Δz₊ₕ`: 深度和厚度数组
- `N`: 土壤层数
- `ibeg`: 模拟起始层索引（`z_bound_top` 对应的下一层）
- `itop`: 上边界层在观测数据中的索引

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

### 6.3 基于配置的高阶 API

位于 `src/Model_Soil.jl`，简化批量模拟 workflow：

| 函数                                                          | 说明                                                                                  |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `setup(config, data_obs)`                                     | 初始化土壤对象和状态，返回 `(soil, state, θ_top, yobs)`                               |
| `Soil_predict(config, theta, state, θ_top)`                   | 运行预测模拟，返回 `ysim`。支持 `method_solve`: `"Bonan"`, `"ODE"`, `"BEPS"`        |
| `Soil_goal(config, theta, state, θ_top, yobs)`                | 计算优化目标函数值（用于 SCE-UA）                                                     |
| `Soil_main(config, data_obs, SITE, dates; maxn, plot_fun...)` | 完整的优化和可视化流程，包括参数优化、结果保存和绘图                                  |

**setup 工作流程**：
1. 调用 `init_SM(config)` 创建土壤对象
2. 使用 `observe_to_state` 将观测数据映射到模型状态
3. 提取上边界层数据 (`θ_top`) 和下层观测数据 (`yobs`)
4. 返回初始化后的对象

**Soil_predict 支持的求解器**：
- `"Bonan"`: 使用 `solve_SM_Bonan` 求解
- `"ODE"`: 使用 `solve_SM_ODE` 求解（需要加载 OrdinaryDiffEqTsit5）
- `"BEPS"`: 使用 `soil_moisture_BEPS` 求解

### 6.4 求解器接口（SoilMoisture/Solve_SM.jl）

| 函数                                                   | 说明                                              |
| ------------------------------------------------------ | ------------------------------------------------- |
| `solve_SM_Bonan(soil, θ_top)`                          | Bonan 隐式求解器，逐时间步求解 Richards 方程      |
| `solve_SM_ODE(soil, θ_top; solver, reltol, abstol)`    | ODE 求解器接口，支持 Tsit5, Rosenbrock23 等       |
| `ModSim_SM(soil, θ_top; method, kw...)`                | 统一求解器接口，method 可选 "Bonan" 或 "ODE"      |
| `SM_param2theta(soil)`                                 | 将 soil 参数转换为优化向量 theta                  |
| `SM_UpdateParam!(soil, theta)`                         | 用优化向量 theta 更新 soil 参数                   |
| `SM_paramBound(soil)`                                  | 返回参数优化上下界 `(lower, upper)`               |

### 6.5 工具函数

位于 `src/ultilize.jl`：

| 函数                             | 说明                   |
| -------------------------------- | ---------------------- |
| `interp_data_depths(A, z, zout)` | 按深度插值观测数据矩阵 |
| `init_grid(zs_center)`           | 根据层中心深度初始化网格 |
| `find_layer_indices(z_bound_top, zs_sim, zs_obs)` | 查找边界层对应的模拟层和观测层索引 |

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

### 7.5 使用 YAML 配置运行（新 API）

```julia
using SoilDifferentialEquations, Ipaper, RTableTools
include("../main_plot.jl")

# 加载配置
fileConfig = "examples/SM/case_SM_China.yaml"
config = load_config(fileConfig)

# 加载并插值观测数据
d = fread(joinpath(dirname(fileConfig), config.file))
data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* config.scale_factor, config.zs_obs_orgin, config.zs_obs)

# 初始化土壤对象和状态
soil, state, θ_top, yobs = setup(config, data_obs)

# 运行完整优化流程（包括可视化）
soil, theta_opt, best_cost = Soil_main(config, data_obs, SITE, dates; 
    plot_fun=plot_result, method_retention="van_Genuchten", maxn=10000)
```

**手动控制优化流程**（可选）：
```julia
# 参数边界和初始值
lower, upper = SM_paramBound(soil)
theta0 = SM_param2theta(soil)

# 手动运行 SCE-UA 优化
theta_opt, feval, _ = sceua(
    theta -> Soil_goal(config, theta, state, θ_top, yobs),
    theta0, lower, upper; maxn=config.maxn)

# 运行预测
ysim = Soil_predict(config, theta_opt, state, θ_top)
```

配置文件示例 (`case_SM_China.yaml`):
```yaml
# SM_China 土壤水模拟配置文件
data:
  file: "data/SM_J1193.csv"     # 数据文件路径（相对于配置文件）
  col_time: 1                   # 时间列索引
  col_obs_start: 2              # 观测数据起始列
  scale_factor: 0.01            # 数据缩放因子

  z_bound_top: 10               # 上边界层深度 [cm]
  zs_obs_orgin: [10, 20, 30, 40, 50, 60, 80, 100]  # 原始观测深度
  zs_obs: [10, 20, 30, 40, 50, 60, 70, 80, 90, 100] # 插值后深度

model:
  soil_type: 7                             # 土壤类型
  same_layer: false                        # 是否各层使用相同参数
  method_retention: "van_Genuchten"        # 土壤持水曲线方法
  method_solve: "Bonan"                    # 求解方法: "Bonan" | "ODE" | "BEPS"
  dt: 3600.0                               # 时间步长（秒）
  zs_center: [2.5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

optimization:
  enable: true      # 是否启用优化
  maxn: 10000       # 最大迭代次数
  objective: "NSE"  # 目标函数: "NSE" | "KGE"

output:
  plot_file: "images/plot_final.png"  # 图表保存路径，留空则不绘图
```

### 7.6 USCRN 站点模拟示例

USCRN (U.S. Climate Reference Network) 数据模拟示例 (`case_SM_Bonan_uscrn.jl`):

```julia
using SoilDifferentialEquations, Ipaper, RTableTools, Dates, LazyArtifacts
include("../main_plot.jl")

fileConfig = "examples/SM/case_SM_Bonan_uscrn.yaml"
config = load_config(fileConfig)

# 加载 USCRN 数据（通过 artifact）
function load_uscrn_data(site_index::Int)
    vars_SM = [:SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, 
               :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]
    f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
    df = fread(f_uscrn2024)
    sites = unique_sort(df.site)
    SITE = sites[site_index]
    d = df[df.site.==SITE, [:time; vars_SM]]
    return d, SITE
end

d, SITE = load_uscrn_data(3)

# 数据准备
(; zs_obs_orgin, zs_obs, scale_factor) = config
data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)

# 使用新 API 运行优化
Soil_main(config, data_obs, SITE, d.time; plot_fun=plot_result, maxn=10000)
```

### 7.7 BEPS 求解器示例

BEPS (Boreal Ecosystem Productivity Simulator) 土壤水分求解器示例，支持自定义土壤参数：

配置文件示例 (`case_SM_BEPS_uscrn.yaml`):
```yaml
# SM_uscrn BEPS 土壤水模拟配置文件

data:
  file: "USCRN_hourly_2024_sp54_Apr-Jun.csv"
  col_time: 1
  col_obs_start: 3       # 土壤水分数据起始列（跳过P_CALC）
  scale_factor: 1.0

  z_bound_top: 5
  zs_obs_orgin: [5, 10, 20, 50, 100]
  zs_obs: [5, 10, 20, 50, 100]

model:
  soil_type: 7
  same_layer: true                       # BEPS 使用统一参数
  method_retention: "van_Genuchten"
  method_solve: "BEPS"                   # 使用 BEPS 求解器
  dt: 3600.0
  zs_center: [1.25, 5, 10, 20, 50, 100]

  # 自定义土壤参数
  soil_params:
    θ_sat: 0.30
    θ_res: 0.03
    Ksat: 1.04
    α: 0.036
    n: 1.56

optimization:
  enable: true
  maxn: 2000
  objective: "KGE"

output:
  plot_file: "plot_BEPS.png"
```

**关键特性：**
- `method_solve: "BEPS"` 启用 BEPS 求解器
- `soil_params` 字段用于覆盖标准土壤类型参数
- BEPS 使用内部迭代直至收敛，通常需要较少的 SCE-UA 迭代次数

### 7.8 土壤温度模拟示例（Tsoil）

```yaml
# Tsoil_CUG 土壤温度模拟配置文件
model:
  model_type: "Tsoil"         # 模型类型: 土壤温度模拟
  soil_type: 7
  dt: 3600.0
  zs_center: [0, 5, 10, 15, 20, 40, 80, 160, 320]
  method_solve: "Bonan"

data:
  file: "../../data/TS_CUG_202306.csv"
  col_time: 1
  col_obs_start: 2
  scale_factor: 1.0
  z_bound_top: 5              # 上边界层深度 (cm)

optimization:
  enable: true
  maxn: 5000
  objective: "NSE"

output:
  plot_file: "plot_Tsoil_CUG.png"
```

### 7.9 土壤质地分类（USDA）

使用 `USDA` 模块进行土壤质地分类：

```julia
using SoilDifferentialEquations.USDA

# 通过砂粒、粉粒含量判断质地类型
texture = soil_texture(60.0, 20.0)  # → SANDY_LOAM (9)

# 枚举类型常量
# CLAY(1), SILTY_CLAY(2), SANDY_CLAY(3), CLAY_LOAM(4), SILTY_CLAY_LOAM(5),
# SANDY_CLAY_LOAM(6), LOAM(7), SILTY_LOAM(8), SANDY_LOAM(9), SILT(10),
# LOAMY_SAND(11), SAND(12)

# 通过名称获取枚举（支持中英文）
parse_soil_type(7)              # → LOAM
parse_soil_type("LOAM")         # → LOAM
parse_soil_type("壤土")          # → LOAM
parse_soil_type(:SILTY_LOAM)    # → SILTY_LOAM

# 通过ID获取枚举
soil_type_from_id(7)  # → LOAM
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
