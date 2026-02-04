# CliMA/ClimaLand.jl 命名规范学习笔记

> 学习来源: https://github.com/CliMA/ClimaLand.jl  
> 分析模块: `src/standalone/Soil/`, `src/shared_utilities/`

---

## 1. 边界命名规范 (Boundary Naming)

### 1.1 核心原则

ClimaLand 明确区分 **Top** 和 **Bottom** 边界，避免模糊的 `boundary` 一词单独使用。

| 命名              | 含义                 | 代码示例                                |
| ----------------- | -------------------- | --------------------------------------- |
| `TopBoundary`     | 顶边界（地表）       | `::ClimaLand.TopBoundary`               |
| `BottomBoundary`  | 底边界（地下深处）   | `::ClimaLand.BottomBoundary`            |
| `TopOfAtmosphere` | 大气顶层（辐射计算） | 辐射模块                                |
| `Surface`         | 地表界面             | `surface_height`, `surface_temperature` |

### 1.2 我们的映射

```yaml
# ❌ 避免
boundary_depth_cm: 10  # 模糊：顶还是底？

# ✅ 推荐（ClimaLand 风格）
z_bound_top: 10       # 地表输入层
top_boundary_depth_cm: 10  # 顶边界深度

# 未来扩展
bot_boundary_depth_cm: 100  # 底边界深度（自由排水/固定水头等）
```

---

## 2. 边界条件类型 (Boundary Condition Types)

### 2.1 State vs Flux 分离

ClimaLand 严格区分**状态边界条件**和**通量边界条件**：

```julia
# 状态边界条件（Dirichlet）
struct MoistureStateBC{F} <: AbstractWaterBC  # θ_l = f(p,t)
struct TemperatureStateBC{F} <: AbstractHeatBC # T = f(p,t)

# 通量边界条件（Neumann）
struct WaterFluxBC{F} <: AbstractWaterBC      # 水分通量
struct HeatFluxBC{F} <: AbstractHeatBC        # 热通量

# 特殊边界条件
struct FreeDrainage <: AbstractWaterBC         # 自由排水
struct AtmosDrivenFluxBC <: AbstractWaterBC    # 大气驱动通量
```

### 2.2 命名模式

| 类型     | 命名模式               | 示例                                    |
| -------- | ---------------------- | --------------------------------------- |
| 状态 BC  | `{Variable}StateBC`    | `MoistureStateBC`, `TemperatureStateBC` |
| 通量 BC  | `{Variable}FluxBC`     | `WaterFluxBC`, `HeatFluxBC`             |
| 抽象类型 | `Abstract{Variable}BC` | `AbstractWaterBC`, `AbstractHeatBC`     |

---

## 3. 空间位置命名 (Spatial Location)

### 3.1 Center vs Face

ClimaLand 使用有限体积法，严格区分**层中心**和**层界面**：

```julia
# 层中心（状态变量位置）
ψ_c        # 中心点水势
temp_c     # 中心点温度
ϑ_l        # 中心点液态水含量

# 层界面（通量计算位置）
ψ_face     # 界面水势
K_face     # 界面导水率
flux_face  # 界面通量

# 空间转换函数
top_center_to_surface(p.soil.ψ)  # 顶层中心 → 地表
```

### 3.2 我们的映射

| 当前命名 | ClimaLand 风格     | 含义             |
| -------- | ------------------ | ---------------- |
| `z`      | `z_centers`        | 层中心深度       |
| `z₊ₕ`    | `z_faces`          | 层界面深度       |
| `Δz`     | `dz`               | 层厚度           |
| `θ`      | `ϑ_l` 或 `theta_l` | 液态水含量       |
| `ψ`      | `ψ_c` / `ψ_face`   | 水势（区分位置） |
| `K`      | `K_face`           | 界面导水率       |

---

## 4. 变量命名风格 (Variable Naming Style)

### 4.1 希腊字母 vs 拉丁字母

ClimaLand 混合使用：

```julia
# 希腊字母（数学/物理标准符号）
ϑ_l        # theta_l: 体积液态水含量
ψ          # psi: 土壤水势
ν          # nu: 孔隙度
θ_r        # theta_r: 残余含水量

# 拉丁字母（描述性变量）
hcm        # hydrology conductivity model
S_s        # specific storage
ν_bc       # boundary condition porosity
```

### 4.2 下标规范

| 下标    | 含义               | 示例           |
| ------- | ------------------ | -------------- |
| `_c`    | center (层中心)    | `ψ_c`          |
| `_face` | face (层界面)      | `K_face`       |
| `_bc`   | boundary condition | `θ_bc`, `ν_bc` |
| `_sfc`  | surface            | `T_sfc`        |
| `_l`    | liquid             | `ϑ_l`, `q_l`   |
| `_liq`  | liquid (verbose)   | `q_liq`        |
| `_ice`  | ice                | `θ_ice`        |

---

## 5. 参数命名 (Parameter Naming)

### 5.1 土壤水力参数

```julia
# 在 ClimaLand 中
struct RichardsParameters{FT, HCM, A, B}
    "基于压力的导水率模型"
    hydrology_cm::HCM
    
    "残余含水量 (m³/m³)"
    θ_r::A  # 或写作 theta_r
    
    "孔隙度/饱和含水量 (m³/m³)"
    ν::A    # nu，等同于 θ_sat
    
    "比储水量 (1/m)"
    S_s::A  # specific storage
    
    "VG 参数: α (1/m)"
    vg_α::A
    
    "VG 参数: n (-)"
    vg_n::A
end
```

### 5.2 我们的映射

| 我们的命名 | ClimaLand 风格 | 说明              |
| ---------- | -------------- | ----------------- |
| `θ_sat`    | `ν`            | 孔隙度/饱和含水量 |
| `θ_res`    | `θ_r`          | 残余含水量        |
| `Ksat`     | `K_sat`        | 饱和导水率        |
| `α`        | `vg_α`         | VG 参数 α         |
| `n`        | `vg_n`         | VG 参数 n         |

---

## 6. 配置文件 YAML 映射

### 6.1 推荐结构 (ClimaLand 风格)

```yaml
# grid.yaml (网格配置)
grid:
  # 垂直坐标定义
  z_centers_cm: [2.5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]  # 层中心
  
  # 边界条件配置
  boundary_conditions:
    top:
      type: "state"           # state | flux | atmosphere_driven
      depth_cm: 10            # 输入数据深度
      variable: "theta_l"     # theta_l | psi | flux
    bottom:
      type: "free_drainage"   # free_drainage | state | flux | impermeable
      # depth_cm: 100         # 可选，底边界深度

  # 模拟范围
  sim_start_layer_idx: 3      # 模拟起始层（默认：top_layer_idx + 1）

# model.yaml (模型配置)
model:
  soil:
    texture: "LOAM"           # 或 "壤土" | 7
    hydrology_model: "van_genuchten"  # van_genuchten | campbell
    
    parameters:  # 可选覆盖
      theta_r: 0.078
      nu: 0.40          # 即 θ_sat
      vg_alpha: 0.036   # 1/cm
      vg_n: 1.56
      K_sat: 1.04       # cm/h

# boundary_data.yaml (边界数据)
boundary_data:
  top:
    file: "SM_J1193.csv"
    time_col: 1
    depth_cm: 10
    variable: "theta_l"       # 变量名
    scale_factor: 0.01
```

### 6.2 对比表

| 当前 (SoilDiffEqs)  | 建议 (ClimaLand 风格)                       | 备注           |
| ------------------- | ------------------------------------------- | -------------- |
| `boundary_depth_cm` | `boundary_conditions.top.depth_cm`          | 嵌套结构更清晰 |
| `method_retention`  | `hydrology_model`                           | 更通用         |
| `method_solve`      | `solver.type`                               | 可扩展         |
| `same_layer`        | `parameters.spatial_variability: "uniform"` | 更描述性       |
| `ibeg`              | `sim_start_layer_idx`                       | 全称           |
| `θ_surf`            | `theta_l_top` 或 `theta_l_surface`          | 明确变量+位置  |

---

## 7. 函数/方法命名 (Function Naming)

### 7.1 ClimaLand 风格

```julia
# 动词 + 名词 + 可选修饰
update_aux!()           # 更新辅助变量
compute_exp_tendency!() # 计算显式趋势
boundary_flux!()        # 计算边界通量
make_update_aux()       # 构造更新函数（工厂模式）

# 空间操作
top_center_to_surface() # 顶层中心 → 地表
```

### 7.2 命名模式总结

| 模式           | 用途                 | 示例                    |
| -------------- | -------------------- | ----------------------- |
| `make_*`       | 工厂函数（构造闭包） | `make_update_aux`       |
| `update_*!`    | 更新操作（mutating） | `update_aux!`           |
| `compute_*!`   | 计算操作             | `compute_exp_tendency!` |
| `get_*`        | 获取器               | `get_domain`            |
| `initialize_*` | 初始化               | `initialize_prognostic` |

---

## 8. 总结：核心命名原则

### 8.1 优先级排序

1. **明确边界**: 使用 `top`/`surface` 和 `bottom`，避免模糊的 `boundary`
2. **区分空间**: `center` vs `face`，`surface` vs `subsurface`
3. **State vs Flux**: `MoistureStateBC` vs `WaterFluxBC`
4. **物理含义**: 使用标准符号 `ϑ_l`, `ψ`, `ν`
5. **下标清晰**: `_c`, `_face`, `_bc`, `_sfc`

### 8.2 快速对照卡

```yaml
# 垂直结构
z_centers: [...]        # 层中心深度
z_faces: [...]          # 层界面深度
dz: [...]               # 层厚度

# 边界条件
boundary_conditions:
  top:
    type: "state"       # state | flux | atmosphere_driven
    depth: 10
    variable: "theta_l"
  bottom:
    type: "free_drainage"

# 土壤参数
theta_r: 0.05           # 残余含水量
nu: 0.45                # 孔隙度 (饱和含水量)
K_sat: 1.0              # 饱和导水率
vg_alpha: 0.036         # VG α
vg_n: 1.56              # VG n

# 模拟设置
sim_start_layer_idx: 3  # 模拟起始层
solver:
  type: "bonan"         # bonan | ode
  dt: 3600.0
```

---

## 参考文件

- `src/standalone/Soil/boundary_conditions.jl` - 边界条件定义
- `src/standalone/Soil/Soil.jl` - 土壤模型核心
- `src/shared_utilities/models.jl` - 通用模型接口
- `src/shared_utilities/boundary_conditions.jl` - 边界类型抽象
