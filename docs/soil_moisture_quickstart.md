# 土壤水运动模拟「一键无脑版」

> 目标：从 **数据准备 → 方案选择 → 参数率定 → 绘图**，用最少步骤跑通。

---

## 0. 你需要准备什么（放到同一文件夹）

1) **观测数据表**（csv/tsv 均可）
- 列格式：
  ```text
  time, P_CALC, SM_5, SM_10, SM_20, SM_50, SM_100
  ```
- 其中 `SM_5` 等是含水量（m3/m3 或体积含水量）

2) **配置文件**（`config_Bonan_CLM5.yaml`）
- 我已经准备好的模板可直接用：
  `test/SM_uscrn/config_Bonan_CLM5.yaml`

---

## 1. 数据准备（最简单）

- 把你的观测文件命名为：`obs.csv`
- 确保观测列从第 3 列开始（time + P_CALC 之后）
- 深度单位：**cm**，如 5,10,20,50,100

配置文件里需要改的只有：
```yaml
model:
  method_retention: ["Campbell", "van_Genuchten"]
  layer_scheme: "two_layer"   # 推荐默认
  two_step: true
  ibeg: 2                      # 默认跳过5cm

optimization:
  enable: true
  maxn: 10000

data:
  file: "obs.csv"
  depths_cm: [5,10,20,50,100]
```

---

## 2. 方案选择（直接抄）

### ✅ 推荐方案
- **Campbell：full 或 two_steps→full**
- **van Genuchten：two_layers 最稳**

只需改一行：
```yaml
layer_scheme: "two_layer"   # 或 full / exp_k / trend
```

---

## 3. 参数率定（无脑跑）

直接运行：
```bash
julia --project test/SM_uscrn/case_solve_Bonan_CLM5.jl
```

输出里会自动显示优化过程与结果：
```text
Optimizing (SCE-UA) Campbell...
Best Cost = -0.78  → NSE ≈ 0.78
```

---

## 4. 绘图（简单版）

最简做法：在 Julia 里跑完后加一段：

```julia
using Plots
plot(obs[:,1], θ_obs[:,1], label="Obs 5cm")
plot!(obs[:,1], ysim[:,1], label="Sim 5cm")
```

> 若需要多深度对比，可循环深度列画多条线。

---

## 5. 一句话流程回顾

✅ **准备观测表 → 修改配置 → 运行脚本 → 看NSE → 绘图对比**

---

## 6. 常见问题（最短版）

**Q1: 结果很差？**
- 先用 Campbell + full / two_layers
- van Genuchten 通常不稳定

**Q2: 运行报错？**
- 检查 obs.csv 路径是否正确
- depths_cm 是否匹配列数

---

如果需要我把这份写成飞书文档/英文版，我可以直接生成。
