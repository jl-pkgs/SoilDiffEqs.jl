#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, size: 13pt, footer: "CUG水文气象学2025", header: "")


// #set page(paper: "a4", margin: 2cm)
// #set text(font: "New Computer Modern", size: 13pt)
#set par(justify: true, leading: 0.65em)
// #set heading(numbering: "1.1")

#show raw: set par(leading: 0.1em, spacing: 0.2em)
#show heading.where(level: 1): set text(size: 18pt, weight: "bold")



#align(center)[
  #text(size: 20pt, weight: "bold")[
    Green-Ampt 下渗模型实现方案
  ]

  // #v(0.5em)
  #text(size: 12pt)[SoilDiffEqs.jl 下渗模块设计]

  // #v(0.3em)
  #text(size: 10pt, style: "italic")[
    日期: 2025-11-04
  ]
]

#v(-1em)

= 1 理论基础

== 1.1 Green-Ampt 模型假设

Green-Ampt (1911) 模型基于以下理想化假设：

1. *均质土壤*：土壤性质在空间上均匀分布
2. *活塞流*：湿润锋为一锐利界面，将土壤分为饱和区和初始含水区
3. *恒定吸力*：湿润锋处保持恒定的基质势 $psi_f$
4. *表层积水*：地表有薄层积水，表层土壤达到饱和

== 1.2 模型推导过程

=== 1.2.1 第一步：建立物理模型

#v(-0.5em)
#figure(
  image("../images/示意图_Green-Ampt.png", width: 60%),
  caption: [
    Green-Ampt模型剖面示意图。\
    定义坐标系$z$向下为负，地表为$z = 0$；湿润锋$z = -L$。湿润锋深度$L$随时间增加。\
    饱和区厚度 $L = F / (Delta theta)$（正值），湿润锋处吸力$psi_f < 0$（负值）。
  ],
) <fig_>

在 Green-Ampt 模型中，土壤剖面被分为两个区域：

- *饱和区* ($-L lt.eq z lt.eq 0$)：含水量为 $theta_s$（饱和）
- *初始区* ($z lt -L$)：含水量为 $theta_i$（初始）



=== 1.2.2 第二步：应用 Darcy 定律

Darcy 定律描述饱和土壤中的水流：

$ q = -K (partial h) / (partial z) $

- $q$：水通量（下渗率）[cm/h]
- $K$：导水率，在饱和区为 $K_s$ [cm/h]
- $h = psi + z$：总水势（基质势 + 重力势）[cm]

=== 1.2.3 第三步：计算水势梯度

*总水势差*：

$ Delta h = h_0 - h_f = 0 - (psi_f + L) = -psi_f - L $


$ (partial h) / (partial z) approx (Delta h) / (Delta z) = (-psi_f - L) / L = -(psi_f + L) / L $

$ f = K_s (1 + psi_f / L) $ <eq:darcy_L>


定义累积下渗量 $F(t)$ [cm]，累积下渗量与湿润锋深度遵循质量守恒。

$ F(t) = integral_0^t f(tau) d tau $


$ F = L dot Delta theta ==> L = F / (Delta theta) $ <eq:L_F>

其中 $Delta theta = theta_s - theta_i$ 是含水量增量。将式 @eq:L_F 代入式 @eq:darcy_L：

$ f = K_s (1 + (psi_f Delta theta) / F) $

这就是 *Green-Ampt 下渗率方程*。

由于 $f = (d F) / (d t)$，得到微分方程：

$ (d F) / (d t) = K_s (1 + (psi_f Delta theta) / F) $ <eq:ga_ode>

将式 @eq:ga_ode 改写为：

$ (d F) / (d t) = K_s + (K_s psi_f Delta theta) / F $

分离变量：

$ F d F = (K_s F + K_s psi_f Delta theta) d t $

定义 $S = psi_f Delta theta$（吸力项），则：

$ (F d F) / (K_s F + K_s S) = d t $

$ (F d F) / (F + S) = K_s d t $ <eq_inf_main>

*对左边进行部分分式分解*：

$ (F) / (F + S) = 1 - S / (F + S) $


#block(
  fill: rgb("#f0f8ff"),
  inset: 10pt,
  radius: 4pt,
  [

    *积分详解*：

    对两边积分：

    $ integral (F d F) / (F + S) = integral (1 - S / (F + S)) d F $

    将右边拆分为两个积分：

    $ = integral 1 d F - integral (S) / (F + S) d F $

    $ = integral 1 d F - S integral (1) / (F + S) d F $

    计算每个积分：
    - 第一项：$integral 1 d F = F + C_1$
    - 第二项：$integral (1)/(F+S) d F = ln(F+S) + C_2$

      合并：

      $ = F + C_1 - S [ln(F + S) + C_2] $

      $ = F - S ln(F + S) + (C_1 - S C_2) $

      将常数合并为一个积分常数 $C = C_1 - S C_2$：

      $ = F - S ln(F + S) + C $

      *关键*：这里的 $C$ 是不定积分产生的积分常数，需要通过初始条件确定。

      #line(length: 100%)

      *类比理解*：

      考虑简单的例子：$integral 2x d x = x^2 + C$

      为什么会有 $C$？因为求导时常数项消失：
      - $(x^2)' = 2x$
      - $(x^2 + 5)' = 2x$
      - $(x^2 - 3)' = 2x$

      所以从 $2x$ 反推原函数时，有无穷多个可能：$x^2 + C$

      在我们的问题中：
      - 不定积分产生：$F - S ln(F+S) + C$
      - 初始条件 $F(0)=0$ 确定：$C = -S ln(S)$
      - 得到唯一解（特解）
  ],
)

因此：

$ F - S ln(F + S) = K_s t + C' $

其中 $C'$ 是右边积分的常数（为了区分，暂记为 $C'$）。

整理后：

$ F - S ln(F + S) = K_s t + C $

=== 1.2.4 第八步：确定积分常数

从上一步得到的通解：

$ F - S ln(F + S) = K_s t + C $

这是含有未知常数 $C$ 的通解。我们需要通过*初始条件*来确定 $C$。

*初始条件*：在下渗开始时刻 $t = 0$，累积下渗量 $F = 0$

将 $t = 0$，$F = 0$ 代入通解：

$ 0 - S ln(0 + S) = K_s dot 0 + C $

$ -S ln(S) = C $

因此：$C = -S ln(S)$

将 $C$ 代回通解：

$ F - S ln(F + S) = K_s t + (-S ln(S)) $

$ F - S ln(F + S) = K_s t - S ln(S) $

#block(
  fill: rgb("#e8ffe8"),
  inset: 10pt,
  radius: 4pt,
  [
    *为什么需要初始条件？*

    微分方程的通解包含积分常数，代表无穷多个可能的解。通过给定初始条件（$t=0$ 时 $F=0$），我们可以从无穷多解中选出唯一的特解，即符合物理实际的那个解。

    *物理意义*：
    - 通解：所有可能的下渗过程（不同的初始状态）
    - 初始条件：从"零时刻开始下渗"这个特定场景
    - 特解：我们真正关心的下渗过程
  ],
)

整理：

$ F + S ln(S / (F + S)) = K_s t $

$ F - S ln((F + S) / S) = K_s t $

$ F - S ln(1 + F / S) = K_s t $

展开 $S = psi_f Delta theta$：

$ F - psi_f Delta theta ln(1 + F / (psi_f Delta theta)) = K_s t $ <eq:ga_integral>

这是 *Green-Ampt 累积下渗量的隐式方程*。

=== 1.2.5 推导总结

从 Darcy 定律出发，通过以下步骤得到 Green-Ampt 方程：

#table(
  columns: 3,
  [*步骤*], [*关键方程*], [*物理意义*],
  [1. Darcy 定律], [$q = -K (partial h) / (partial z)$], [饱和流基本定律],
  [2. 水势梯度], [$(partial h) / (partial z) = -(psi_f + L) / L$], [表层到湿润锋的势差],
  [3. 下渗率], [$f = K_s (1 + psi_f / L)$], [与湿润锋深度相关],
  [4. 质量守恒], [$F = L Delta theta$], [累积下渗 = 湿润层蓄水],
  [5. 最终形式], [$f = K_s (1 + (psi_f Delta theta) / F)$], [Green-Ampt 方程],
  [6. 积分解], [式 @eq:ga_integral], [时间-累积下渗关系],
)

#block(
  fill: rgb("#e8f4f8"),
  inset: 10pt,
  radius: 4pt,
  [
    *推导核心要点*

    Green-Ampt 模型的推导建立在三个关键假设上：

    1. *活塞流假设*：允许将复杂的水势分布简化为阶跃函数
      - 饱和区：$theta = theta_s$，$K = K_s$
      - 初始区：$theta = theta_i$

    2. *平均水势梯度*：将整个饱和区的水势梯度近似为平均值
      - 避免求解复杂的 Richards 方程
      - 用 $Delta h / L$ 代替 $partial h / partial z$

    3. *质量守恒*：建立累积下渗量 $F$ 与湿润锋深度 $L$ 的关系
      - $F = L dot Delta theta$
      - 将空间变量转化为时间变量

    *推导的巧妙之处*：通过质量守恒关系 $L = F / Delta theta$，将含有空间变量 $L$ 的 Darcy 方程转化为只含时间积分量 $F$ 的常微分方程，从而避免了求解偏微分方程。

    *数学简化*：
    $ f = K_s (1 + psi_f / L) quad arrow.r.double quad f = K_s (1 + (psi_f Delta theta) / F) $

    这一转化是 Green-Ampt 模型的核心创新。
  ],
)

== 1.3 控制方程

基于 Darcy 定律，下渗率 $f$ (infiltration rate) 为：

$ f = K_s [1 + (psi_f Delta theta) / F] $

其中：
- $K_s$：饱和导水率 [cm/h]
- $psi_f$：湿润锋吸力水头 [cm]（正值）
- $Delta theta = theta_s - theta_i$：含水量差
- $F$：累积下渗量 [cm]

累积下渗量随时间演化：

$ (d F) / (d t) = f = K_s [1 + (psi_f Delta theta) / F] $

== 1.4 物理意义解释

=== 1.4.1 下渗率的两个组成部分

Green-Ampt 方程可以改写为：

$ f = K_s + (K_s psi_f Delta theta) / F $

$ f = f_"gravity" + f_"suction" $

其中：
- *重力项* $f_"gravity" = K_s$：由重力驱动的稳态下渗
- *吸力项* $f_"suction" = (K_s psi_f Delta theta) / F$：由基质势梯度驱动的额外下渗

=== 1.4.2 下渗率的时间演化

观察方程可知：

#table(
  columns: 2,
  [*时间阶段*], [*下渗特征*],
  [$t -> 0$，$F -> 0$], [$f -> infinity$（吸力项占主导）],
  [$t$ 增大，$F$ 增大], [$f$ 递减（吸力项减小）],
  [$t -> infinity$，$F -> infinity$], [$f -> K_s$（重力项占主导）],
)

物理解释：

1. *初期*：湿润锋浅（$L$ 小），水势梯度大，下渗快
2. *中期*：湿润锋加深（$L$ 增大），梯度减小，下渗减慢
3. *后期*：吸力影响可忽略，以重力排水为主

=== 1.4.3 与湿润锋深度的关系

由 $L = F / (Delta theta)$ 可得：

$ f = K_s (1 + psi_f / L) $

说明：
- 湿润锋越浅（$L$ 小），下渗率越大
- 湿润锋越深（$L$ 大），下渗率越接近 $K_s$

=== 1.4.4 参数的物理意义

#table(
  columns: 3,
  [*参数*], [*物理含义*], [*对下渗的影响*],
  [$K_s$ 大], [土壤渗透性好], [下渗快，长期稳态下渗率高],
  [$psi_f$ 大], [湿润锋吸力强], [初期下渗快，但随时间衰减],
  [$Delta theta$ 大], [土壤初始很干], [蓄水能力强，下渗快],
  [$theta_i$ 小], [初始含水量低], [增大 $Delta theta$，下渗增强],
)

== 1.5 积水时间 $t_p$

=== 1.5.1 物理过程

当降雨强度 $P > K_s$ 时，下渗过程分为两个阶段：

*阶段 1：积水前* ($0 < t < t_p$)
- 地表尚未形成积水
- 表层土壤未达到饱和
- 所有降雨都能下渗：$f = P$

*阶段 2：积水后* ($t > t_p$)
- 地表开始积水
- 表层土壤饱和
- 下渗率由 Green-Ampt 方程控制：$f = K_s (1 + (psi_f Delta theta) / F)$

积水发生的临界条件：下渗能力等于降雨率。

=== 1.5.2 积水时累积下渗量 $F_p$ 推导

在积水时刻 $t_p$，下渗率从降雨控制转为土壤控制，满足：

$ f(t_p) = P $

根据 Green-Ampt 方程，在积水时刻：

$ f = K_s (1 + (psi_f Delta theta) / F_p) = P $

求解 $F_p$：

$ K_s + (K_s psi_f Delta theta) / F_p = P $

$ F_p = (K_s psi_f Delta theta) / (P - K_s) $ <eq:Fp>

#block(
  fill: rgb("#fff0f5"),
  inset: 10pt,
  radius: 4pt,
  [
    *物理意义*：

    - 分子 $K_s psi_f Delta theta$：土壤吸力产生的"额外下渗潜力"
    - 分母 $P - K_s$：超出饱和导水率的"多余降雨"
    - $F_p$ 越大：需要更多时间才积水（土壤吸力强或多余降雨少）

    *特殊情况*：
    - 当 $P -> K_s$，$F_p -> infinity$（永不积水）
    - 当 $P >> K_s$，$F_p -> (K_s psi_f Delta theta) / P$（很快积水）
  ]
)

=== 1.5.3 积水时间 $t_p$ 推导

在积水前阶段，下渗率恒定为 $f = P$，因此：

$ F_p = integral_0^(t_p) f d t = integral_0^(t_p) P d t = P dot t_p $

求解 $t_p$：

$ t_p = F_p / P $ <eq:tp_basic>

将式 @eq:Fp 代入式 @eq:tp_basic：

$ t_p = (K_s psi_f Delta theta) / (P(P - K_s)) $ <eq:tp>

#block(
  fill: rgb("#e8ffe8"),
  inset: 10pt,
  radius: 4pt,
  [
    *推导步骤总结*：

    1. 积水条件：$K_s(1 + psi_f Delta theta / F_p) = P$

    2. 求解 $F_p$：$F_p = (K_s psi_f Delta theta) / (P - K_s)$

    3. 积水前质量守恒：$F_p = P dot t_p$

    4. 得到：$t_p = F_p / P = (K_s psi_f Delta theta) / (P(P - K_s))$
  ]
)

=== 1.5.4 积水时间的影响因素

从式 @eq:tp 可以分析各参数的影响：

#table(
  columns: 3,
  [*参数增大*], [*对 $t_p$ 的影响*], [*物理解释*],
  [$K_s$ 增大], [$t_p$ 增大], [土壤渗透性好，不易积水],
  [$psi_f$ 增大], [$t_p$ 增大], [吸力强，延迟积水],
  [$Delta theta$ 增大], [$t_p$ 增大], [土壤干燥，蓄水能力强],
  [$P$ 增大], [$t_p$ 减小], [降雨强度大，快速积水],
)

*数值例子*：

假设 $K_s = 1.0$ cm/h，$psi_f = 10.0$ cm，$Delta theta = 0.3$

#table(
  columns: 3,
  [*P (cm/h)*], [*$F_p$ (cm)*], [*$t_p$ (h)*],
  [1.5], [6.0], [4.0],
  [2.0], [3.0], [1.5],
  [3.0], [1.5], [0.5],
  [5.0], [0.75], [0.15],
)

观察：降雨强度越大，积水越快

== 1.6 积水后的下渗

由式#[@eq_inf_main]，$S = psi_f Delta theta$，从积水时刻 $t_p$（累积下渗 $F_p$）到时刻 $t$（累积下渗 $F$）积分：

$ integral_(F_p)^F (F d F) / (F + S) = K_s integral_(t_p)^t d t $

*左边积分*（利用 1.2.3 节的结果）：

$ integral_(F_p)^F (F d F) / (F + S) = [F - S ln(F + S)]_(F_p)^F $

$ = [F - S ln(F + S)] - [F_p - S ln(F_p + S)] $

$ = (F - F_p) - S [ln(F + S) - ln(F_p + S)] $

$ = (F - F_p) - S ln((F + S) / (F_p + S)) $

*右边积分*：

$ K_s integral_(t_p)^t d t = K_s (t - t_p) $

因此得到隐式方程：

$ (F - F_p) - S ln((F + S) / (F_p + S)) = K_s (t - t_p) $

展开 $S = psi_f Delta theta$：

$ boxed(K_s (t - t_p) = (F - F_p) - psi_f Delta theta ln[(F + psi_f Delta theta) / (F_p + psi_f Delta theta)]) $ <eq:ga_post_ponding>

#block(
  fill: rgb("#e8f4f8"),
  inset: 10pt,
  radius: 4pt,
  [
    *与完整方程的关系*：

    完整 Green-Ampt 方程（从 $t=0, F=0$ 积分）：
    $ K_s t = F - psi_f Delta theta ln(1 + F / (psi_f Delta theta)) $

    积水后方程（从 $t_p, F_p$ 积分）：
    $ K_s (t - t_p) = (F - F_p) - psi_f Delta theta ln[(F + psi_f Delta theta) / (F_p + psi_f Delta theta)] $

    两者的区别在于初始条件不同。
  ]
)

#pagebreak()

= 2 数值实现方法

== 2.1 方法 1：显式时间步进（推荐用于小时间步长）

*算法流程*：

```julia
# 每个时间步 dt：
if !ponding_occurred
    # 阶段 1：积水前
    if P <= Ksat
        f = P  # 下渗率等于降雨率
        F += f * dt
    else
        # 检查是否达到积水
        if F >= F_p
            ponding_occurred = true
            ponding_time = t
        else
            f = P
            F += f * dt
        end
    end
else
    # 阶段 2：积水后
    f = Ksat * (1 + ψf * Δθ / F)
    F += f * dt
end
```

*优点*：
- 简单直观
- 易于实现
- 适合小时间步长（dt < 1分钟）

*缺点*：
- 大时间步长可能不稳定

== 2.2 方法 2：隐式求解（适用于大时间步长）

使用 Newton-Raphson 迭代求解隐式方程。

对于时间步 $n -> n+1$：

$ K_s Delta t = F^(n+1) - F^n - psi_f Delta theta ln[(psi_f Delta theta + F^(n+1)) / (psi_f Delta theta + F^n)] $

定义残差函数：

$
  R(F^(n+1)) = F^(n+1) - F^n - K_s Delta t - psi_f Delta theta ln[(psi_f Delta theta + F^(n+1)) / (psi_f Delta theta + F^n)]
$

导数：

$ (d R) / (d F^(n+1)) = 1 - (psi_f Delta theta) / (psi_f Delta theta + F^(n+1)) $

Newton 迭代：

$ F^(n+1)_(k+1) = F^(n+1)_k - R(F^(n+1)_k) / ((d R) / (d F^(n+1)))|_(F^(n+1)_k) $

*收敛判据*：$|R(F^(n+1))| < epsilon$，取 $epsilon = 10^(-6)$ cm

== 2.3 方法 3：显式-隐式混合

```julia
# 显式预测
f_pred = Ksat * (1 + ψf * Δθ / F_n)
F_pred = F_n + f_pred * dt

# 隐式校正
f_corr = Ksat * (1 + ψf * Δθ / F_pred)
F_new = F_n + 0.5 * (f_pred + f_corr) * dt
```

= 3 Julia 实现代码

#block(
  fill: rgb("#e8f4f8"),
  inset: 10pt,
  radius: 4pt,
  width: 15cm,
  [
    *代码文件位置*

    完整实现代码已分离到独立文件，便于维护和使用：

    - 核心算法：`src/Infiltration/GreenAmpt.jl`

    - 参数估算：`src/Infiltration/Parameters.jl`

    - 混合模型：`src/Infiltration/Hybrid.jl`

    - 主模块：`src/Infiltration/Infiltration.jl`

    - 测试：`test/test-infiltration.jl`

    - 示例: `examples/example_greenampt.jl`

    以下仅展示关键代码片段，完整代码请查看上述文件。
  ],
)

= 4 与 Richards 方程耦合

== 4.1 混合模型策略

在土壤剖面模拟中，可以采用混合策略：
- *表层*：使用 Green-Ampt（计算效率高）
- *深层*：使用 Richards 方程（精度高）

```julia
"""
混合下渗模型
"""
mutable struct HybridInfiltration{FT}
    ga::GreenAmpt{FT}           # Green-Ampt 模型
    soil::Soil{FT}              # Richards 方程土壤
    coupling_depth::FT          # 耦合深度 [m]
    use_GA::Bool                # 是否使用 GA
end

function update_hybrid_infiltration!(
    hi::HybridInfiltration,
    P::FT,
    dt::FT
) where FT

    if hi.use_GA
        # 使用 Green-Ampt
        f = infiltration_GreenAmpt!(hi.ga, P, dt)

        # 将下渗率作为 Richards 方程顶部边界
        hi.soil.Q0 = -f

        # 检查表层是否达到饱和
        if hi.ga.F > hi.coupling_depth * 100  # 转换为 cm
            # 切换到完整 Richards 方程
            hi.use_GA = false
            hi.soil.ψ0 = 0.0  # 饱和边界
        end

    else
        # 使用完整 Richards 方程
        cal_Q!(hi.soil; Q0=-P, method="Q0")
    end

    return nothing
end
```

= 5 参数估算方法

== 5.1 从土壤质地估算

使用 Rawls et al. (1983) 的经验关系：

#table(
  columns: 5,
  [*质地*], [$K_s$ (cm/h)], [$theta_s$], [$theta_i$], [$psi_f$ (cm)],
  [砂土], [11.78], [0.437], [0.078], [4.95],
  [壤质砂土], [2.99], [0.437], [0.090], [6.13],
  [砂壤土], [2.18], [0.453], [0.125], [11.01],
  [壤土], [1.32], [0.463], [0.155], [8.89],
  [粉壤土], [0.68], [0.501], [0.210], [16.68],
  [砂质黏壤土], [0.43], [0.398], [0.186], [21.85],
  [黏壤土], [0.23], [0.464], [0.242], [20.88],
  [粉质黏壤土], [0.15], [0.471], [0.255], [27.30],
  [黏土], [0.06], [0.475], [0.286], [31.63],
)

== 5.2 从 Van Genuchten 参数转换

```julia
function estimate_GA_parameters(vg::ParamVanGenuchten{FT}) where FT
    (; θ_sat, θ_res, Ksat, α, n, m) = vg

    # 方法 1: Morel-Seytoux (1978)
    ψf_1 = 0.5 / α

    # 方法 2: Bouwer (1966) - 使用进气值
    ψf_2 = 1.0 / α

    # 方法 3: 数值积分（更精确）
    # ψf = ∫[θi to θs] ψ(θ) dθ / Δθ
    θ_range = range(θ_res + 0.01, θ_sat - 0.01, length=100)
    ψ_values = Retention_ψ.(θ_range; par=vg)
    ψf_3 = -trapz(θ_range, ψ_values) / (θ_sat - θ_res)

    return (ψf_MS = ψf_1, ψf_B = ψf_2, ψf_exact = ψf_3)
end
```

= 6 测试与验证

== 6.1 单元测试

```julia
@testset "GreenAmpt Model" begin
    # 测试 1: 小降雨，无积水
    ga = GreenAmpt(Ksat=1.0, θi=0.1, θs=0.4, ψf=10.0)
    P = 0.5  # cm/h
    dt = 600.0  # 10 分钟

    f = infiltration_GreenAmpt!(ga, P, dt)
    @test f ≈ P  # 下渗率应等于降雨率
    @test !ga.ponding_occurred

    # 测试 2: 大降雨，积水
    ga2 = GreenAmpt(Ksat=1.0, θi=0.1, θs=0.4, ψf=10.0)
    P2 = 5.0  # cm/h，超过 Ksat

    for i = 1:100
        f = infiltration_GreenAmpt!(ga2, P2, dt, method=:implicit)
    end

    @test ga2.ponding_occurred
    @test ga2.F > 0
    @test ga2.f < P2  # 积水后下渗率小于降雨率

    # 测试 3: 质量守恒
    ga3 = GreenAmpt(Ksat=2.0, θi=0.15, θs=0.45, ψf=15.0)
    P3 = 3.0
    dt3 = 60.0
    ntim = 1000

    total_infiltration = 0.0
    for i = 1:ntim
        f = infiltration_GreenAmpt!(ga3, P3, dt3)
        total_infiltration += f * dt3/3600
    end

    @test abs(total_infiltration - ga3.F) < 1e-6
end
```

== 6.2 与解析解对比

对于恒定降雨 $P > K_s$，积水后的理论解为：

$ F(t) = K_s t + psi_f Delta theta ln[1 + F(t)/(psi_f Delta theta)] $

```julia
function analytical_solution_GA(Ksat, ψf, Δθ, t)
    # Newton 求解
    F = Ksat * t
    for iter = 1:20
        R = F - Ksat * t - ψf * Δθ * log(1 + F/(ψf*Δθ))
        dR = 1 - ψf * Δθ / (ψf*Δθ + F)
        F -= R / dR
    end
    return F
end

@testset "Compare with analytical solution" begin
    ga = GreenAmpt(Ksat=1.5, θi=0.1, θs=0.4, ψf=12.0)
    P = 10.0  # 大降雨，立即积水
    dt = 60.0
    times = 0.0:0.1:10.0  # 10 小时

    # 数值解
    F_numerical = Float64[]
    for t in times
        if t > 0
            infiltration_GreenAmpt!(ga, P, dt, method=:implicit)
        end
        push!(F_numerical, ga.F)
    end

    # 解析解
    F_analytical = analytical_solution_GA.(ga.Ksat, ga.ψf, ga.Δθ, times)

    # 对比
    max_error = maximum(abs.(F_numerical .- F_analytical))
    @test max_error < 0.1  # cm
end
```

= 7 性能优化

== 7.1 时间步长控制

```julia
function adaptive_timestep_GA(ga::GreenAmpt, P, dt_max)
    """
    自适应时间步长选择
    """
    (; f, F, Ksat, ψf, Δθ) = ga

    # 计算局部截断误差估计
    df_dF = -Ksat * ψf * Δθ / F^2

    # CFL 条件
    dt_cfl = 0.5 / abs(df_dF)

    # 取最小值
    dt = min(dt_max, dt_cfl)

    return max(dt, 1.0)  # 至少 1 秒
end
```

== 7.2 向量化处理

对于多列土壤剖面：

```julia
struct GreenAmptField{FT}
    models::Vector{GreenAmpt{FT}}
    ncols::Int
end

function infiltration_GreenAmpt_vectorized!(
    gaf::GreenAmptField,
    P::Vector,
    dt
)
    f_out = similar(P)

    Threads.@threads for i = 1:gaf.ncols
        f_out[i] = infiltration_GreenAmpt!(
            gaf.models[i],
            P[i],
            dt
        )
    end

    return f_out
end
```

= 8 总结与建议

== 8.1 Green-Ampt 优势

1. *计算效率*：比 Richards 方程快 100-1000 倍
2. *参数少*：只需 4 个参数
3. *物理意义清晰*：易于理解和调试
4. *稳定性好*：不存在 Richards 方程的数值困难

== 8.2 Green-Ampt 局限性

1. *假设严格*：要求均质土壤、活塞流
2. *初始条件*：要求均匀初始含水量
3. *再分布*：无法模拟下渗后的水分再分布
4. *分层土壤*：难以处理多层土壤

== 8.3 应用建议

*使用 Green-Ampt 的场景*：
- 短时强降雨事件模拟
- 大尺度水文模型（需要高效率）
- 实时洪水预报
- 参数敏感性分析

*使用 Richards 方程的场景*：
- 精确的土壤水分动态
- 长时间模拟（包括干湿循环）
- 分层/非均质土壤
- 根系吸水耦合

*混合策略*：
- 降雨期间用 Green-Ampt
- 降雨后用 Richards 方程模拟再分布

== 8.4 下一步工作

1. 实现基础 Green-Ampt 模块
2. 编写完整测试套件
3. 与 Richards 方程对比验证
4. 性能基准测试
5. 文档和示例

= 9 参考文献

+ Green, W.H., Ampt, G.A. (1911). Studies on soil physics, part 1: The flow of air and water through soils. _Journal of Agricultural Science_, 4(1), 1-24.

+ Rawls, W.J., Brakensiek, D.L., Miller, N. (1983). Green-Ampt infiltration parameters from soils data. _Journal of Hydraulic Engineering_, 109(1), 62-70.

+ Morel-Seytoux, H.J., Khanji, J. (1974). Derivation of an equation of infiltration. _Water Resources Research_, 10(4), 795-800.

+ Bouwer, H. (1966). Rapid field measurement of air entry value and hydraulic conductivity of soil as significant parameters in flow system analysis. _Water Resources Research_, 2(4), 729-738.

+ Mein, R.G., Larson, C.L. (1973). Modeling infiltration during a steady rain. _Water Resources Research_, 9(2), 384-394.
