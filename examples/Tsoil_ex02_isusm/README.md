# ISUSM 土壤温度模拟案例 (Case 02)

本案例展示了如何使用 `SoilDiffEqs.jl` 模拟土壤温度剖面，并结合 ISUSM 站点观测数据进行参数优化（率定）。

提供了两种运行方式：
1.  **配置文件模式 (推荐)**: 通过修改 `config.yaml` 即可运行，无需编写代码。
2.  **脚本模式 (高级)**: 直接修改 Julia 脚本，适合深度定制。

---

## 方式一：配置文件模式 (推荐)

### 1. 环境准备
使用此模式需要安装 `YAML.jl` 包。请在 Julia REPL 中运行：
```julia
using Pkg
Pkg.add("YAML")
```

### 2. 快速运行
直接运行驱动脚本即可：

```bash
# 默认使用同目录下的 config.yaml
julia --project examples/Tsoil_ex02_isusm/run_config.jl

# 或指定配置文件路径
julia --project examples/Tsoil_ex02_isusm/run_config.jl examples/Tsoil_ex02_isusm/config.yaml
```

### 3. 配置说明 (`config.yaml`)
你可以复制 `examples/Tsoil_ex02_isusm/config.yaml` 并根据需要修改：

```yaml
# 数据配置
data:
  file: "data/isusm_TS_202207.csv" # 输入CSV路径
  time_col: "time"                 # 时间列名
  obs_start_col: 4                 # 观测数据起始列索引
  time_steps: 240                  # 模拟时长（时间步数）

# 模拟设置
simulation:
  dt: 3600.0                       # 时间步长(s)

# 模型网格
model:
  dz: 0.05                         # 网格间距(m)
  depths_inch: [4, 12, ... 52]     # 观测深度(英寸), 脚本会自动转为米
  Tsurf_layer_index: 3             # 选取第几层观测作为上边界(Tsurf)
  soil_type: 7                     # 土壤类型参数

# 优化设置
optimization:
  enable: true                     # 是否开启参数率定
  method: "SCE-UA"
  max_iterations: 50000            # 最大迭代次数
  bounds:                          # 参数范围
    kappa: [0.1, 10.0]
    cv: [1.0e6, 5.0e6]

# 输出
output:
  plot: true
  save_fig: "result.png"
```

---

## 方式二：脚本模式 (高级)

直接编辑并运行 `case02_Tsoil_isusm.jl`。

### 1. 准备数据
脚本需要读取 CSV 格式的观测数据文件。默认路径为 `$dir_soil/data/isusm_TS_202207.csv`。

**数据格式要求**:
*   格式: CSV (逗号分隔)
*   必需列:
    *   `time`: 时间列
    *   第4列及之后: 各层土壤温度观测值

### 2. 运行模型
在 Julia REPL 中运行以下命令：

```julia
# 激活项目环境
using Pkg; Pkg.activate(".")

# 运行脚本
include("examples/Tsoil_ex02_isusm/case02_Tsoil_isusm.jl")
```

## 模型参数说明

*   **优化变量**:
    *   `κ` (kappa): 土壤热导率 (Thermal Conductivity)。
    *   `cv`: 土壤体积热容量 (Volumetric Heat Capacity)。
*   **目标函数**:
    *   最大化观测值与模拟值的 Nash-Sutcliffe 效率系数 (NSE)。
