"""
# 地表积水处理方案测试示例

本文件演示如何使用三种地表积水处理方案：
- 方案A：基于蓄水容量
- 方案B：基于坡面汇流（曼宁公式）
- 方案C：混合方案

## 使用方法

```julia
# 在Julia REPL中运行
include("dev/ponding_schemeA.jl")
include("dev/ponding_schemeB.jl")
include("dev/ponding_schemeC.jl")
include("dev/test_ponding_schemes.jl")

# 运行测试
test_all_schemes()
```
"""

# 加载三个方案
include("ponding_schemeA.jl")
include("ponding_schemeB.jl")
include("ponding_schemeC.jl")


"""
    test_schemeA()

测试方案A：蓄水容量方案
"""
function test_schemeA()
  println("\n" * "="^60)
  println("测试方案A：基于蓄水容量的方案")
  println("="^60)

  # 场景：平坦农田
  h_pond_max = 0.5  # [cm]
  dt = 0.1  # [h]

  println("\n场景：平坦农田 (h_pond_max = $h_pond_max cm)")
  println("-"^60)

  # 测试1：小雨，不产流
  println("\n测试1：小雨 (P=2 cm/h, inf=1.5 cm/h)")
  h_current = 0.0
  P, inf = 2.0, 1.5
  runoff, h_new = ponding_schemeA(h_current, P, inf, dt; h_pond_max=h_pond_max)
  println("  初始积水: $(h_current) cm")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")

  # 测试2：中雨，开始产流
  println("\n测试2：中雨 (P=8 cm/h, inf=3 cm/h)")
  h_current = 0.3
  P, inf = 8.0, 3.0
  runoff, h_new = ponding_schemeA(h_current, P, inf, dt; h_pond_max=h_pond_max)
  println("  初始积水: $(h_current) cm")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")

  # 测试3：暴雨，大量产流
  println("\n测试3：暴雨 (P=30 cm/h, inf=2 cm/h)")
  h_current = 0.2
  P, inf = 30.0, 2.0
  runoff, h_new = ponding_schemeA(h_current, P, inf, dt; h_pond_max=h_pond_max)
  println("  初始积水: $(h_current) cm")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")
end


"""
    test_schemeB()

测试方案B：坡面汇流方案
"""
function test_schemeB()
  println("\n" * "="^60)
  println("测试方案B：基于坡面汇流的方案（曼宁公式）")
  println("="^60)

  # 场景：坡地农田
  S0 = 0.05  # 5% 坡度
  n_manning = 0.15
  L_char = 50.0  # [m]
  h_pond_min = 0.1  # [cm]
  dt = 0.1  # [h]

  println("\n场景：坡地农田 (S0=$S0, n=$n_manning, L=$L_char m)")
  println("-"^60)

  # 测试1：积水未达到起始产流深度
  println("\n测试1：积水较少 (h < h_pond_min)")
  h_current = 0.05
  P, inf = 5.0, 3.0
  runoff, h_new = ponding_schemeB(h_current, P, inf, dt;
                                  S0=S0, n_manning=n_manning, L_char=L_char, h_pond_min=h_pond_min)
  println("  初始积水: $(h_current) cm")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")

  # 测试2：中等积水深度
  println("\n测试2：中等积水 (h=1.0 cm)")
  h_current = 1.0
  P, inf = 5.0, 2.0
  runoff, h_new = ponding_schemeB(h_current, P, inf, dt;
                                  S0=S0, n_manning=n_manning, L_char=L_char, h_pond_min=h_pond_min)
  # 计算排水速率
  h_excess = h_current - h_pond_min
  h_m = h_excess / 100.0
  v = (1.0 / n_manning) * (h_m^(2/3)) * sqrt(S0)
  q = v * h_m * 3600.0 * 100.0
  drain_rate = q / L_char
  println("  初始积水: $(h_current) cm")
  println("  计算排水速率: $(round(drain_rate, digits=2)) cm/h")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")

  # 测试3：较大积水深度
  println("\n测试3：较大积水 (h=2.0 cm)")
  h_current = 2.0
  P, inf = 5.0, 2.0
  runoff, h_new = ponding_schemeB(h_current, P, inf, dt;
                                  S0=S0, n_manning=n_manning, L_char=L_char, h_pond_min=h_pond_min)
  # 计算排水速率
  h_excess = h_current - h_pond_min
  h_m = h_excess / 100.0
  v = (1.0 / n_manning) * (h_m^(2/3)) * sqrt(S0)
  q = v * h_m * 3600.0 * 100.0
  drain_rate = q / L_char
  println("  初始积水: $(h_current) cm")
  println("  计算排水速率: $(round(drain_rate, digits=2)) cm/h")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")
  println("  注意：h从1.0→2.0 cm，排水速率增加 $(round(drain_rate/5.39, digits=2)) 倍 (h^(5/3)效应)")
end


"""
    test_schemeC()

测试方案C：混合方案
"""
function test_schemeC()
  println("\n" * "="^60)
  println("测试方案C：混合方案（三阶段）")
  println("="^60)

  # 场景：通用场景
  h_pond_min = 0.1  # [cm]
  h_pond_max = 2.0  # [cm]
  S0 = 0.01  # 1% 坡度
  n_manning = 0.15
  L_char = 100.0  # [m]
  dt = 0.1  # [h]

  println("\n场景：通用场景")
  println("  h_pond_min=$h_pond_min cm, h_pond_max=$h_pond_max cm")
  println("  S0=$S0, n=$n_manning, L=$L_char m")
  println("-"^60)

  # 测试1：阶段1 - 微地形蓄水
  println("\n测试1：阶段1 - 微地形蓄水 (h < h_pond_min)")
  h_current = 0.05
  P, inf = 3.0, 2.5
  runoff, h_new = ponding_schemeC(h_current, P, inf, dt;
                                  h_pond_min=h_pond_min, h_pond_max=h_pond_max,
                                  S0=S0, n_manning=n_manning, L_char=L_char)
  println("  初始积水: $(h_current) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")
  println("  状态：微地形蓄水，无径流")

  # 测试2：阶段2 - 坡面汇流
  println("\n测试2：阶段2 - 坡面汇流 (h_pond_min < h < h_pond_max)")
  h_current = 1.0
  P, inf = 5.0, 2.0
  runoff, h_new = ponding_schemeC(h_current, P, inf, dt;
                                  h_pond_min=h_pond_min, h_pond_max=h_pond_max,
                                  S0=S0, n_manning=n_manning, L_char=L_char)
  println("  初始积水: $(h_current) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")
  println("  状态：坡面汇流")

  # 测试3：阶段3 - 超蓄产流
  println("\n测试3：阶段3 - 超蓄产流 (h > h_pond_max)")
  h_current = 1.5
  P, inf = 20.0, 2.0
  runoff, h_new = ponding_schemeC(h_current, P, inf, dt;
                                  h_pond_min=h_pond_min, h_pond_max=h_pond_max,
                                  S0=S0, n_manning=n_manning, L_char=L_char)
  println("  初始积水: $(h_current) cm")
  println("  净增水量: $((P-inf)*dt) cm")
  println("  新积水深度: $(round(h_new, digits=3)) cm")
  println("  径流: $(round(runoff, digits=3)) cm")
  println("  状态：超蓄产流")
end


"""
    compare_schemes()

比较三种方案在相同条件下的表现
"""
function compare_schemes()
  println("\n" * "="^60)
  println("比较三种方案")
  println("="^60)

  # 相同的初始条件
  h_current = 0.5  # [cm]
  P = 10.0  # [cm/h]
  inf = 3.0  # [cm/h]
  dt = 0.1  # [h]

  println("\n初始条件：")
  println("  积水深度: $h_current cm")
  println("  降水强度: $P cm/h")
  println("  下渗速率: $inf cm/h")
  println("  时间步长: $dt h")
  println("-"^60)

  # 方案A
  runoff_A, h_new_A = ponding_schemeA(h_current, P, inf, dt; h_pond_max=0.5)
  println("\n方案A (h_pond_max=0.5 cm):")
  println("  新积水: $(round(h_new_A, digits=3)) cm")
  println("  径流: $(round(runoff_A, digits=3)) cm")

  # 方案B
  runoff_B, h_new_B = ponding_schemeB(h_current, P, inf, dt;
                                      S0=0.05, n_manning=0.15, L_char=50.0, h_pond_min=0.1)
  println("\n方案B (S0=0.05, n=0.15, L=50m):")
  println("  新积水: $(round(h_new_B, digits=3)) cm")
  println("  径流: $(round(runoff_B, digits=3)) cm")

  # 方案C
  runoff_C, h_new_C = ponding_schemeC(h_current, P, inf, dt;
                                      h_pond_min=0.1, h_pond_max=2.0,
                                      S0=0.01, n_manning=0.15, L_char=100.0)
  println("\n方案C (h_min=0.1, h_max=2.0, S0=0.01, L=100m):")
  println("  新积水: $(round(h_new_C, digits=3)) cm")
  println("  径流: $(round(runoff_C, digits=3)) cm")

  println("\n对比分析：")
  println("  方案A产流最多（蓄水容量小，立即超蓄）")
  println("  方案B产流中等（坡度大，排水较快）")
  println("  方案C产流最少（容量大+坡度小，调蓄能力强）")
end


"""
    simulate_rainfall_event()

模拟一次完整的降雨事件
"""
function simulate_rainfall_event()
  println("\n" * "="^60)
  println("模拟降雨事件：间歇性降雨")
  println("="^60)

  # 降雨过程：30分钟暴雨 → 30分钟停雨 → 30分钟中雨
  dt = 1/6  # 10分钟一个时间步 [h]
  P_series = [20.0, 20.0, 20.0,  # 前30分钟：暴雨
              0.0, 0.0, 0.0,     # 中间30分钟：停雨
              8.0, 8.0, 8.0]     # 后30分钟：中雨
  inf_series = [5.0, 4.0, 3.5,   # 下渗速率逐渐降低
                2.0, 2.0, 2.0,
                3.0, 2.5, 2.0]

  println("\n使用方案C模拟：")
  println("时间步长: $(dt*60) 分钟")
  println("-"^60)
  println("时刻(min)  降雨(cm/h)  下渗(cm/h)  积水(cm)  径流(cm)")
  println("-"^60)

  h_pond = 0.0
  total_runoff = 0.0

  for i in 1:length(P_series)
    P = P_series[i]
    inf = inf_series[i]

    runoff, h_pond = ponding_schemeC(h_pond, P, inf, dt;
                                     h_pond_min=0.1, h_pond_max=2.0,
                                     S0=0.02, n_manning=0.15, L_char=50.0)
    total_runoff += runoff

    t_min = i * dt * 60
    println("$(lpad(round(Int, t_min),8))   $(lpad(P,10))   $(lpad(inf,10))   " *
            "$(lpad(round(h_pond, digits=3),8))  $(lpad(round(runoff, digits=3),8))")
  end

  println("-"^60)
  println("总径流: $(round(total_runoff, digits=3)) cm")
  println("总降雨: $(round(sum(P_series) * dt, digits=3)) cm")
  println("径流系数: $(round(total_runoff / (sum(P_series) * dt), digits=3))")
end


"""
    test_all_schemes()

运行所有测试
"""
function test_all_schemes()
  println("\n╔" * "="^58 * "╗")
  println("║  地表积水处理方案测试程序" * " "^30 * "║")
  println("╚" * "="^58 * "╝")

  test_schemeA()
  test_schemeB()
  test_schemeC()
  compare_schemes()
  simulate_rainfall_event()

  println("\n" * "="^60)
  println("所有测试完成！")
  println("="^60)
  println("\n提示：")
  println("  1. 方案A适合平坦地区，参数简单")
  println("  2. 方案B适合坡地，考虑地形影响")
  println("  3. 方案C通用性强，推荐使用")
  println("  4. 参数需根据具体场景率定")
  println()
end


# 如果直接运行此文件，执行所有测试
if abspath(PROGRAM_FILE) == @__FILE__
  test_all_schemes()
end
