#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#let Qin = "Q_in"
#let Qout = "Q_out"


== 方案1：下渗模型

$ (Qin - Qout) Delta t = Delta theta_(t,1) Delta z_1 $

对第一层土壤进行建模：

- Q_in 代表流入第一层土壤的水量，下渗量；

- Q_out 代表流出第一层土壤的水量。

$
  Qout = - K pdv(psi + z, z)
  = - K_1 ([psi_2 + z_2] - [psi_1 + z_1]) / (Delta z_2)
$

#box-red()[
  - $K_1$: 第一层土壤的水力传导系数未知。
  - $psi$: 土壤水力参数，5~6个参数。
]

== 方案2：整层方案


