#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#import "@preview/tablex:0.0.9": cellx, colspanx, rowspanx, tablex

// 代码块样式
#show raw.where(block: true): it => block(
  width: 100%,
  fill: rgb("#f5f5f5"),
  inset: 10pt,
  radius: 4pt,
  text(size: 9.5pt, it),
)

// 标题页
#align(center)[
  #v(3cm)
  #text(size: 28pt, weight: "bold")[
    SoilDiffEqs.jl
  ]

  #v(0.5cm)
  #text(size: 20pt)[
    模型使用手册
  ]

  #v(1cm)
  #text(size: 14pt, style: "italic")[
    土壤水热运动微分方程求解包
  ]

  #v(3cm)
  #text(size: 12pt)[
    版本：v1.0 \
    更新日期：2025-12-19
  ]
]

#pagebreak()

// 目录
#outline(
  title: [目录],
  indent: auto,
  depth: 3,
)

#pagebreak()

= 1 模型概述

*SoilDiffEqs.jl* 是一个用于求解土壤水热运动微分方程的Julia软件包。该模型能够模拟土壤中的水分运动、热量传输以及地下水动态，为生态水文、农业灌溉、气候模拟等领域提供重要的计算工具。

== 1.1 主要功能

- *土壤水分模拟*：基于Richards方程求解土壤水分的三维运动
- *土壤温度模拟*：基于热传导方程计算土壤温度剖面
- *地下水模拟*：模拟地下水位动态及其与土壤水分的相互作用
- *植被-土壤耦合*：考虑根系吸水和蒸散发过程

== 1.2 模型特点

- *高效求解*：提供多种数值求解方案，计算速度快
- *灵活的边界条件*：支持通量边界和水势边界两种边界条件
- *多种土壤参数化方案*：支持van Genuchten和Campbell两种土壤水力特性模型
- *模块化设计*：各功能模块独立，便于扩展和维护

#pagebreak()

= 2 物理背景与理论基础

== 2.1 土壤水分运动

=== 2.1.1 Richards方程的推导

土壤水分运动由*Richards方程*描述，该方程结合了达西定律和质量守恒原理：

$ frac(partial theta, partial t) = frac(partial, partial z) [K(theta) (frac(partial psi, partial z) + 1)] - S $

其中：
- $theta$ — 体积含水量 [$"m"^3"/"m^3$]
- $t$ — 时间 [s]
- $z$ — 深度坐标（向下为负）[m]
- $K(theta)$ — 水力传导度 [cm/h]
- $psi$ — 土壤水势 [cm]
- $S$ — 源汇项（如根系吸水）[cm/h]

#block(
  fill: rgb("#e8f4f8"),
  inset: 10pt,
  radius: 4pt,
)[
  *物理意义*：
  - *左侧项* $frac(partial theta, partial t)$：表示土壤含水量的局部时间变化率
  - *右侧第一项*：表示水分通量 $q$ 的垂直梯度
    - 水势梯度驱动：$frac(partial psi, partial z)$ 水从高水势流向低水势
    - 重力驱动：常数项 $1$ 代表重力作用，水向下渗透
  - *右侧第二项* $S$：源汇项，如根系吸水（负值）或灌溉（正值）
]

=== 2.1.2 达西定律

水分通量 $q$ 遵循达西定律的扩展形式：

$ q = -K(theta) (frac(partial psi, partial z) + 1) $

*物理解释*：
- 土壤水分流动的驱动力包括两部分：
  + *水势梯度*：毛细力和吸附力引起的扩散
  + *重力*：恒定向下的拉力
- 导水率 $K(theta)$ 随含水量非线性变化：
  + 饱和时：$K = K_"sat"$（最大值）
  + 干燥时：$K arrow 0$（几乎不透水）

=== 2.1.3 水分运动的物理过程

*入渗过程*（降雨后）：
+ *初期*：地表水势高，水势梯度大，入渗快
+ *湿润锋向下推进*：上层逐渐饱和，水势梯度减小
+ *稳定入渗*：当表层接近饱和，入渗率趋于稳定值（接近 $K_"sat"$）

*再分布过程*（降雨停止后）：
+ 重力作用占主导，水分继续向下运动
+ 上层土壤开始脱水，水势下降
+ 水分在剖面上重新分布，趋向平衡状态

*蒸发过程*：
+ 表层水分首先蒸发，表层含水量急剧下降
+ 水势梯度反转，深层水分向上运动
+ 形成"干燥锋"，逐渐向下推进

== 2.2 土壤水力特性

土壤的水力特性描述了土壤含水量、水势和导水率之间的非线性关系。模型支持两种参数化方案：

=== 2.2.1 van Genuchten模型（推荐）

该模型通过一组经验公式描述土壤水力特性：

*水分特征曲线*（水势-含水量关系）：

$ theta(psi) = theta_r + frac(theta_s - theta_r, [1 + |alpha psi|^n]^m) $

*导水率函数*：

$ K(theta) = K_"sat" sqrt(S_e) [1 - (1 - S_e^(1\/m))^m]^2 $

其中：
- $S_e = frac(theta - theta_r, theta_s - theta_r)$ — 有效饱和度
- $theta_s$ — 饱和含水量 [$"m"^3"/"m^3$]
- $theta_r$ — 残余含水量 [$"m"^3"/"m^3$]
- $K_"sat"$ — 饱和导水率 [cm/h]
- $alpha$ — 进气值的倒数 [$"cm"^(-1)$]
- $n, m$ — 孔隙分布指数，通常 $m = 1 - 1\/n$

#block(
  fill: rgb("#fff4e6"),
  inset: 10pt,
  radius: 4pt,
)[
  *参数的物理意义*：

  *饱和含水量* $theta_s$：
  - 所有孔隙被水充满时的含水量
  - 砂土 < 壤土 < 粘土（孔隙度增大）

  *残余含水量* $theta_r$：
  - 强吸附于土粒表面，无法流动的水分
  - 代表土壤持水的"死区"

  *进气参数* $alpha$：
  - 表征土壤开始脱水的难易程度
  - $alpha$ 大：大孔隙多，容易失水（砂土）
  - $alpha$ 小：小孔隙多，持水能力强（粘土）

  *孔隙分布* $n$：
  - 反映孔隙大小的均匀性
  - $n$ 大：孔隙大小差异大，水分特征曲线陡峭
  - $n$ 小：孔隙较均匀，水分特征曲线平缓
]

*曲线形态特征*：
- *饱和区* ($psi arrow 0$)：$theta arrow theta_s$，所有孔隙充满水
- *过渡区* ($-100 < psi < 0$)：含水量快速下降，大孔隙排水
- *干燥区* ($psi < -100$)：含水量缓慢下降，仅小孔隙保持水分
- *残余区* ($psi arrow -infinity$)：$theta arrow theta_r$，仅强结合水

=== 2.2.2 Campbell模型

这是一个更简单的经验模型：

$ theta(psi) = theta_s (frac(psi, psi_"sat"))^(-1\/b) $

$ K(theta) = K_"sat" (frac(theta, theta_s))^(2b+3) $

参数 $b$ 表征土壤孔隙分布，砂土约为4，粘土约为11。

== 2.3 土壤温度传输

=== 2.3.1 热传导方程

土壤温度变化由*热传导方程*控制：

$ C_v frac(partial T, partial t) = frac(partial, partial z) [kappa frac(partial T, partial z)] $

其中：
- $T$ — 土壤温度 [°C]
- $C_v$ — 体积热容 [J/($"m"^3 dot.op$K)]
- $kappa$ — 热导率 [W/(m·K)]

#block(
  fill: rgb("#e8f4f8"),
  inset: 10pt,
  radius: 4pt,
)[
  *物理意义*：
  - *左侧项*：单位体积土壤的热量变化率 = $C_v dot Delta T / Delta t$
  - *右侧项*：热通量梯度，热量从高温流向低温
  - *扩散性质*：与水分运动类似，都是扩散过程
]

=== 2.3.2 土壤热力学参数

*热导率* $kappa$ [W/(m·K)]：
- *定义*：单位温度梯度下的热通量密度
- *物理意义*：衡量土壤传递热量的能力
- *影响因素*：
  + 土壤质地：砂土 > 壤土 > 粘土
  + 含水量：湿土 >> 干土（水 $kappa approx 0.6$，空气 $kappa approx 0.025$）
  + 密实度：压实土壤热导率高

*体积热容* $C_v$ [J/($"m"^3 dot.op$K)]：
- *定义*：单位体积土壤升温1K所需热量
- *物理意义*：衡量土壤储存热量的能力
- *计算公式*：
  $ C_v = theta_"water" rho_"water" c_"water" + theta_"solid" rho_"solid" c_"solid" + theta_"air" rho_"air" c_"air" $
- 其中 $theta$ 为各相体积分数，水的热容远大于空气

=== 2.3.3 温度传输的物理过程

*日变化规律*：
+ *表层*（0-5cm）：
  - 温度波动剧烈，振幅可达15-20°C
  - 白天快速升温，夜间快速降温
  - 响应时间：分钟级

+ *浅层*（5-30cm）：
  - 振幅衰减至5-10°C
  - 相位滞后2-4小时
  - 温度峰值出现在下午

+ *深层*（30-100cm）：
  - 振幅继续衰减至1-3°C
  - 相位滞后4-8小时
  - 温度变化平缓

*年变化规律*：
+ 表层年振幅20-30°C
+ 深层（2-3m）年振幅 $<$ 5°C
+ 深度约10-15m处温度恒定（年平均温度）

*热扩散系数*：
定义为 $D_T = kappa / C_v$ [m²/s]，表征温度扩散速度：
- 砂土：$D_T approx 5 times 10^(-7)$ m²/s
- 壤土：$D_T approx 3 times 10^(-7)$ m²/s
- 粘土：$D_T approx 2 times 10^(-7)$ m²/s

== 2.4 地下水动态

模型采用简化的地下水模型，通过*水位深度 (zwt)* 跟踪地下水位：

$ S_y frac(d z_"wt", d t) = "Recharge" - "Drainage" - "Baseflow" $

其中：
- $z_"wt"$ — 地下水位深度 [m]
- $S_y$ — 给水度 [$"m"^3"/"m^3$]
- Recharge — 入渗补给 [mm/d]
- Drainage — 深层渗漏 [mm/d]
- Baseflow — 基流 [mm/d]

#block(
  fill: rgb("#e8f4f8"),
  inset: 10pt,
  radius: 4pt,
)[
  *物理意义*：地下水位的升降取决于补给和排泄的平衡。当降雨充沛时，入渗补给大于排泄，水位上升；干旱时则相反。
]

#pagebreak()

= 3 土壤水分模块

== 3.1 核心功能

土壤水分模块求解Richards方程，计算土壤各层的含水量随时间的变化。

*主要函数*：
- `soil_moisture!()` — 求解基于水势边界条件的土壤水分
- `soil_moisture_Q0!()` — 求解基于通量边界条件的土壤水分
- `RichardsEquation()` — Richards方程的微分方程形式

== 3.2 边界条件

模型支持两种上边界条件：

=== 3.2.1 第一类边界（水势边界，ψ0）

指定地表水势，通常用于模拟积水或饱和条件：

```julia
ψ0 = 0.0  # 饱和条件，单位：cm
```

=== 3.2.2 第二类边界（通量边界，Q0）

指定入渗速率，常用于降雨-下渗过程：

```julia
Q0 = -2.0  # 入渗速率，负值表示向下，单位：cm/h
```

#block(
  fill: rgb("#e6f3ff"),
  inset: 8pt,
  radius: 4pt,
)[
  *选择建议*：
  - 无积水条件：使用通量边界 `Q0`
  - 有积水或灌溉：使用水势边界 `ψ0 = 0`
]

== 3.3 求解方法

模型提供两种求解方案：

+ *隐式有限差分法*（推荐）
  - 计算稳定
  - 允许较大时间步长
  - 效率高

+ *常微分方程求解器（ODE）*
  - 公式清晰
  - 适合研究用途
  - 计算较慢

== 3.4 使用示例

```julia
using SoilDifferentialEquations

# 创建土壤对象
Δz = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]  # 各层厚度，单位：m
soil = Soil(Δz)

# 设置初始条件
soil.θ .= 0.2  # 初始含水量
soil.dt = 3600.0  # 时间步长，单位：秒

# 设置边界条件（降雨入渗）
Q0 = -1.5  # 入渗速率 1.5 cm/h

# 设置根系吸水（蒸散发）
soil.sink .= 0.01  # 单位：cm/h

# 求解一个时间步
soil_moisture_Q0!(soil, Q0)

# 查看结果
println("土壤含水量：", soil.θ)
println("土壤水势：", soil.ψ, " cm")
```

#pagebreak()

= 4 土壤温度模块

== 4.1 核心功能

土壤温度模块求解热传导方程，模拟土壤温度剖面的时空变化。

*主要函数*：
- `soil_temperature!()` — 已知地表温度求解土壤温度
- `soil_temperature_F0!()` — 已知地表热通量求解土壤温度

== 4.2 边界条件

=== 4.2.1 指定地表温度 (Tsurf)

```julia
Tsurf = 25.0  # 地表温度，单位：°C
Tsoil_next, G = soil_temperature!(soil, Tsurf)
```

=== 4.2.2 指定地表热通量 (F0)

```julia
F0 = -50.0  # 向下的热通量，单位：W/m²
Tsoil_next, Tsurf = soil_temperature_F0!(soil, F0)
```

== 4.3 数值方法

提供两种数值格式：

+ *全隐式格式*（默认，推荐）
  - 无条件稳定
  - 精度略低

+ *Crank-Nicolson格式*
  - 精度更高（二阶）
  - 计算量稍大

*使用方法*：

```julia
# 隐式格式
soil_temperature!(soil, Tsurf; solution="implicit")

# Crank-Nicolson格式
soil_temperature!(soil, Tsurf; solution="crank-nicolson")
```

== 4.4 使用示例

```julia
using SoilDifferentialEquations

# 创建土壤对象
Δz = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
soil = Soil(Δz)

# 设置初始温度
soil.Tsoil .= 20.0  # 初始温度 20°C

# 设置热力参数
soil.param.κ .= 1.5    # 热导率 W/(m·K)
soil.param.cv .= 2.5e6  # 体积热容 J/(m³·K)

# 设置时间步长
soil.dt = 3600.0  # 1小时

# 指定地表温度
Tsurf_next = 30.0  # 地表温度升至30°C

# 求解
Tsoil_next, G = soil_temperature!(soil, Tsurf_next)

# 查看结果
println("土壤温度剖面：", Tsoil_next)
println("地表热通量 G = ", G, " W/m²")
```

#pagebreak()

= 5 地下水模块

== 5.1 功能说明

地下水模块模拟浅层地下水位的动态变化，以及地下水与非饱和土壤水分的相互作用。

*核心概念*：
- *zwt (water table depth)*：地下水位深度，向下为正 [m]
- *jwt (water table index)*：地下水层的索引
- *wa (aquifer storage)*：潜水含水层储水量 [mm]
- *Specific yield (Sy)*：给水度，表征含水层蓄水能力

== 5.2 主要函数

- `find_jwt()` — 根据水位深度确定地下水所在层
- `GW_Update_ZWT!()` — 更新地下水位
- `GW_Correctθ!()` — 修正受地下水影响层的含水量

== 5.3 物理过程

=== 5.3.1 地下水位更新

地下水位根据补给和排泄的平衡更新：

$ Delta z_"wt" = frac("Recharge" - "Drainage", S_y) $

=== 5.3.2 土壤含水量修正

当地下水位上升/下降时，需要修正地下水面以下土壤层的含水量：
- 水位以下：$theta = theta_"sat"$（饱和）
- 水位以上：根据水势平衡计算

== 5.4 使用示例

```julia
using SoilDifferentialEquations

# 创建土壤对象
soil = Soil(Δz)

# 设置初始地下水位
soil.zwt = -2.0  # 地下水位在地表下2米

# 计算地下水层索引
soil.jwt = find_jwt(soil.z₊ₕ, soil.zwt)

# 地下水补给（来自底层渗漏）
recharge = 5.0  # mm/day
dt_day = 1.0    # 1天

# 更新地下水位
GW_Update_ZWT!(soil, recharge, dt_day)

# 修正土壤含水量
GW_Correctθ!(soil)

println("新的地下水位：", soil.zwt, " m")
println("地下水层索引：", soil.jwt)
```

#pagebreak()

= 6 模型参数说明

== 6.1 土壤水力参数

=== 6.1.1 van Genuchten参数（推荐使用）

#figure(
  tablex(
    columns: 7,
    align: center + horizon,
    auto-vlines: false,
    auto-hlines: false,

    // 表头
    rowspanx(2)[*参数*],
    rowspanx(2)[*符号*],
    rowspanx(2)[*单位*],
    colspanx(3)[*典型值*],
    (),
    (),
    (),
    (),
    (),
    [*砂土*],
    [*壤土*],
    [*粘土*],

    [饱和含水量],
    [$theta_"sat"$],
    [$"m"^3"/"m^3$],
    [0.43],
    [0.40],
    [0.38],
    [残余含水量],
    [$theta_"res"$],
    [$"m"^3"/"m^3$],
    [0.045],
    [0.078],
    [0.068],
    [饱和导水率],
    [$K_"sat"$],
    [cm/h],
    [29.7],
    [1.04],
    [0.2],
    [进气参数],
    [$alpha$],
    [$"cm"^(-1)$],
    [0.145],
    [0.036],
    [0.008],
    [孔隙参数],
    [$n$],
    [—],
    [2.68],
    [1.56],
    [1.09],
  ),
  caption: [van Genuchten模型参数典型值],
)

#text(size: 9pt, style: "italic")[
  数据来源：Carsel & Parrish (1988) 土壤分类数据库
]

=== 6.1.2 Campbell参数

#figure(
  tablex(
    columns: 4,
    align: center + horizon,
    auto-vlines: false,

    [*参数*],
    [*符号*],
    [*单位*],
    [*典型范围*],
    [饱和含水量],
    [$theta_"sat"$],
    [$"m"^3"/"m^3$],
    [0.3 - 0.5],
    [饱和水势],
    [$psi_"sat"$],
    [cm],
    [-5 to -30],
    [饱和导水率],
    [$K_"sat"$],
    [cm/h],
    [0.1 - 50],
    [孔隙分布],
    [$b$],
    [—],
    [4 - 13],
  ),
  caption: [Campbell模型参数典型值],
)

== 6.2 土壤热力参数

#figure(
  tablex(
    columns: 5,
    align: center + horizon,
    auto-vlines: false,

    [*参数*],
    [*符号*],
    [*单位*],
    [*干土*],
    [*湿土*],
    [热导率],
    [$kappa$],
    [W/(m·K)],
    [0.2 - 0.5],
    [1.0 - 2.5],
    [体积热容],
    [$C_v$],
    [J/($"m"^3 dot.op$K)],
    [$1.3 times 10^6$],
    [$2.5 times 10^6$],
  ),
  caption: [土壤热力参数范围],
)

#block(
  fill: rgb("#fff4e6"),
  inset: 8pt,
  radius: 4pt,
)[
  *说明*：土壤热力参数受含水量影响显著，湿土热导率和热容均高于干土。
]

== 6.3 植被参数

#figure(
  tablex(
    columns: 4,
    align: center + horizon,
    auto-vlines: false,

    [*参数*],
    [*符号*],
    [*单位*],
    [*说明*],
    [根系分布],
    [$f_"root"$],
    [—],
    [各层根系比例，和为1],
    [蒸散发],
    [ET],
    [cm/h],
    [总蒸散速率],
    [土壤水分胁迫],
    [$beta$],
    [—],
    [0-1之间，越小越干旱],
  ),
  caption: [植被参数说明],
)

#pagebreak()

= 7 使用指南

== 7.1 安装

```julia
using Pkg
Pkg.add(url="https://github.com/jl-pkgs/SoilDifferentialEquations.jl")
```

== 7.2 快速开始

=== 7.2.1 示例1：模拟降雨入渗过程

```julia
using SoilDifferentialEquations

# 1. 设置土壤分层
Δz = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]  # 7层土壤
soil = Soil(Δz; method_retention="van_Genuchten")

# 2. 设置土壤参数（壤土）
soil.param.θ_sat .= 0.40
soil.param.θ_res .= 0.078
soil.param.Ksat .= 1.04
soil.param.α .= 0.036
soil.param.n .= 1.56

# 3. 初始条件
soil.θ .= 0.15  # 初始较干

# 4. 降雨事件
rainfall = 2.0  # cm/h，持续3小时
duration = 3 * 3600  # 秒
soil.dt = 600.0  # 10分钟步长

# 5. 时间循环
nsteps = Int(duration / soil.dt)
for i in 1:nsteps
    Q0 = -rainfall
    soil_moisture_Q0!(soil, Q0)

    # 每小时输出一次
    if mod(i, 6) == 0
        println("时间 $(i*soil.dt/3600) h,
                 平均含水量 = $(mean(soil.θ[1:soil.N]))")
    end
end

println("最终含水量剖面：", soil.θ[1:soil.N])
```

=== 7.2.2 示例2：土壤温度日变化

```julia
using SoilDifferentialEquations

# 1. 创建土壤
Δz = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5]
soil = Soil(Δz)

# 2. 设置热力参数
soil.param.κ .= 1.5
soil.param.cv .= 2.0e6

# 3. 初始条件
soil.Tsoil .= 20.0
soil.dt = 1800.0  # 30分钟

# 4. 模拟24小时
for hour in 1:24
    # 地表温度日变化（正弦曲线）
    Tsurf = 20.0 + 10.0 * sin(2π * hour / 24)

    # 求解2次（1小时）
    for _ in 1:2
        soil_temperature!(soil, Tsurf)
        soil.Tsoil .= soil.u[1:soil.N]  # 更新温度
    end

    println("$(hour):00 - 地表 $(round(Tsurf, digits=1))°C,
             10cm $(round(soil.Tsoil[3], digits=1))°C")
end
```

=== 7.2.3 示例3：考虑地下水的长期模拟

```julia
using SoilDifferentialEquations

# 创建土壤
soil = Soil(Δz; method_retention="van_Genuchten")

# 初始地下水位
soil.zwt = -1.5  # 1.5米深
soil.jwt = find_jwt(soil.z₊ₕ, soil.zwt)

# 初始化平衡含水量
Equilibrium!(soil)

# 长期模拟（30天）
days = 30
for day in 1:days
    # 日降雨（随机）
    P = rand() < 0.3 ? rand(0.5:0.1:3.0) : 0.0

    # 日蒸散发
    ET = 0.3  # cm/day
    soil.sink .= ET / 24.0  # cm/h

    # 每小时求解
    for hour in 1:24
        Q0 = -P
        soil_moisture_Q0!(soil, Q0)

        # 更新地下水
        recharge = soil.Q[end] * soil.dt / 3600 / 10  # mm/day
        GW_Update_ZWT!(soil, recharge, 1/24)
        GW_Correctθ!(soil)
    end

    println("Day $day: 降雨=$(P) cm,
             地下水位=$(round(soil.zwt, digits=2)) m")
end
```

== 7.3 高级用法

=== 7.3.1 自定义土壤参数

```julia
# 使用USDA土壤分类
using SoilDifferentialEquations.USDA

# 获取"sandy loam"的参数
par = USDA.get_soil_param("sandy loam")
soil.param.θ_sat .= par.θ_sat
soil.param.Ksat .= par.Ksat
# ... 其他参数
```

=== 7.3.2 水热耦合模拟

```julia
# 同时求解土壤水分和温度
for step in 1:nsteps
    # 1. 更新土壤水分
    soil_moisture_Q0!(soil, Q0)

    # 2. 根据含水量更新热力参数
    update_thermal_properties!(soil)

    # 3. 求解土壤温度
    soil_temperature!(soil, Tsurf)

    # 4. 更新状态
    soil.θ_prev .= soil.θ
    soil.Tsoil .= soil.u[1:soil.N]
end
```

#pagebreak()

= 8 常见问题

== 8.1 Q1: 如何选择时间步长？

#block(
  fill: rgb("#e6f3ff"),
  inset: 10pt,
  radius: 4pt,
)[
  *建议*：
  - 土壤水分：300 - 3600秒（5分钟 - 1小时）
  - 土壤温度：1800 - 7200秒（30分钟 - 2小时）
  - 降雨入渗：更小步长（60 - 600秒）

  *原则*：步长越小越精确，但计算时间越长。
]

== 8.2 Q2: 模型不收敛怎么办？

*可能原因*：
+ 时间步长过大 → 减小 `dt`
+ 土壤参数不合理 → 检查参数范围
+ 边界条件突变 → 使用渐变边界条件

== 8.3 Q3: 如何验证模型结果？

*建议*：
+ 质量守恒检验：`error_SM(soil)` 检查水量平衡
+ 能量守恒检验：温度求解自动检查
+ 与观测数据对比

== 8.4 Q4: 模型适用范围？

#grid(
  columns: 2,
  gutter: 1em,
  block(
    fill: rgb("#e8f8e8"),
    inset: 10pt,
    radius: 4pt,
  )[
    *适用*：
    - 均质/分层土壤
    - 一维垂直水热传输
    - 中等含水量范围
  ],
  block(
    fill: rgb("#ffe8e8"),
    inset: 10pt,
    radius: 4pt,
  )[
    *不适用*：
    - 优先流（macropore flow）
    - 强冻土相变
    - 二维/三维问题
  ],
)

#pagebreak()

= 9 参考文献

+ *Richards方程*：
  Richards, L.A. (1931). Capillary conduction of liquids through porous mediums. _Physics_, 1(5), 318-333.

+ *van Genuchten模型*：
  van Genuchten, M.T. (1980). A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. _Soil Science Society of America Journal_, 44(5), 892-898.

+ *数值方法*：
  Celia, M.A., et al. (1990). A general mass-conservative numerical solution for the unsaturated flow equation. _Water Resources Research_, 26(7), 1483-1496.

+ *土壤热传导*：
  Bonan, G.B. (2019). _Climate Change and Terrestrial Ecosystem Modeling_. Cambridge University Press.

+ *参数数据库*：
  Carsel, R.F., & Parrish, R.S. (1988). Developing joint probability distributions of soil water retention characteristics. _Water Resources Research_, 24(5), 755-769.
