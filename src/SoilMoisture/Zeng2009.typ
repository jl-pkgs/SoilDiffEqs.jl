#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#counter(heading).update(6)
#set par(leading: 1.24em)
#let delta(x) = $Delta #x$

== 1 CLM土壤水运动求解方案



=== 1.1 $q$的偏微分

$ q = -K pdv(psi + z, z) $ <q>

$
  q_(i-1) & = -K_(i-1/2) ( (psi_i + z_i) - (psi_(i-1) + z_(i-1))) / (z_i - z_(i-1))
$

$
  q_(i) & = -K_(i+1/2) ( (psi_(i+1) + z_(i+1)) - (psi_(i) + z_(i))) / (z_(i+1) - z_(i))
$


对$q$求偏导：

$
  pdv(q_(i-1), theta_(i-1)) =
  K_(i-1/2) / (z_i - z_(i-1)) pdv(psi_(i-1), theta_(i-1))
  - pdv(K_(i-1/2), theta_(i-1)) [(psi_i + z_i - psi_(i-1) - z_(i-1)) / (z_i - z_(i-1))]
$ <dq1>

$
  pdv(q_(i-1), theta_(i)) =
  -K_(i-1/2) / (z_i - z_(i-1)) pdv(psi_(i), theta_(i))
  - pdv(K_(i-1/2), theta_(i)) [(psi_i + z_i - psi_(i-1) - z_(i-1)) / (z_i - z_(i-1))]
$

$
  pdv(q_(i), theta_(i)) =
  K_(i+1/2) / (z_(i+1) - z_(i)) pdv(psi_(i), theta_(i))
  - pdv(K_(i+1/2), theta_(i)) [(psi_(i+1) + z_(i+1) - psi_(i) - z_(i)) / (z_(i+1) - z_(i))]
$

$
  pdv(q_(i), theta_(i+1)) =
  -K_(i+1/2) / (z_(i+1) - z_(i)) pdv(psi_(i+1), theta_(i+1))
  - pdv(K_(i+1/2), theta_(i+1)) [(psi_(i+1) + z_(i+1) - psi_(i) - z_(i)) / (z_(i+1) - z_(i))]
$ <dq4>


=== 1.2 均衡水势$psi_E$

为了考虑土壤水与地下水的交互作用，Zeng et al. (2009)引入了均衡水势的概念$psi_E$。$q = 0$时，对应的土壤水势，即均衡水势$psi_E$，对应的土壤含水量为均衡含水量$theta_E$。当处于均衡状态时，每一层土壤的总势能均相等：

$ psi_E + z = psi_"sat" [theta_E(z) / theta_"sat"]^(-B) + z = psi_"sat" + z_"wt" = C $

$
  q & = -K pdv((psi + z), z) = -K pdv((psi + z - C), z) \
    & = -K pdv((psi - psi_E), z)
$ <q_psi_e>

将公式#[@q_psi_e]应用于离散形式：

$
  q_(i-1) & = -K_(i-1/2) ( (psi_i - psi_(E,i)) - (psi_(i-1) - psi_(E,i-1))) / (z_i - z_(i-1))
$

$
  q_(i) & = -K_(i+1/2) ( (psi_(i+1) - psi_(E,i+1)) - (psi_(i) - psi_(E,i))) / (z_(i+1) - z_(i))
$

由于均衡水势$psi_E$由地下水位$z_"wt"$和土壤性质决定，不随实际含水量$theta$变化，即$pdv(psi_E, theta) = 0$。因此，对$q$求偏导时：

$
  pdv(q_(i-1), theta_(i-1)) =
  K_(i-1/2) / (z_i - z_(i-1)) pdv(psi_(i-1), theta_(i-1))
  - pdv(K_(i-1/2), theta_(i-1)) [(psi_i - psi_(E,i) - psi_(i-1) + psi_(E,i-1)) / (z_i - z_(i-1))]
$

$
  pdv(q_(i-1), theta_(i)) =
  -K_(i-1/2) / (z_i - z_(i-1)) pdv(psi_(i), theta_(i))
  - pdv(K_(i-1/2), theta_(i)) [(psi_i - psi_(E,i) - psi_(i-1) + psi_(E,i-1)) / (z_i - z_(i-1))]
$

$
  pdv(q_(i), theta_(i)) =
  K_(i+1/2) / (z_(i+1) - z_(i)) pdv(psi_(i), theta_(i))
  - pdv(K_(i+1/2), theta_(i)) [(psi_(i+1) - psi_(E,i+1) - psi_(i) + psi_(E,i)) / (z_(i+1) - z_(i))]
$

$
  pdv(q_(i), theta_(i+1)) =
  -K_(i+1/2) / (z_(i+1) - z_(i)) pdv(psi_(i+1), theta_(i+1))
  - pdv(K_(i+1/2), theta_(i+1)) [(psi_(i+1) - psi_(E,i+1) - psi_(i) + psi_(E,i)) / (z_(i+1) - z_(i))]
$



=== 1.3 构造三角矩阵

在微积分中，若$z = f(x, y)$，则$d z = pdv(z, x)d y + pdv(z, y) d x$。若$q_(i-1) = f(theta_(i-1), theta_i)$，根据积分法则可以得到：

$ d q_(i-1) = pdv(q_(i-1), theta_(i-1)) d theta_(i-1) + pdv(q_(i-1), theta_(i)) d theta_(i) $

// q = $#(q)_#iq$
#let dqdw(q, w) = pdv(q, w) + $Delta$ + w

#let calqRight(i0, i1, t0, t1) = {
  $q_#i0^#t0 + dqdw(q_#i0, theta_(#i0)) + dqdw(q_#i0, theta_(#i1))$
}

#let calq(i0, i1, t0, t1) = {
  $q_#i0^#t1 = #calqRight(i0, i1, t0, t1)$
}

$ calq(i-1, i, n, n+1) $ <q_i0>
$ calq(i, i+1, n, n+1) $ <q_i1>

$
  (Delta theta_(i)) / (Delta t) Delta z_i = -q_(i-1)^(n+1) + q_(i)^(n+1) - e_i
$ <d_theta>

其中，上标$n$和$n+1$代表时间；下标$i-1$、$i$、$i+1$代表土壤层。

#beamer-block(
  [
    1. $q$为矢量，向下为负，$-q$是为了从矢量转为标量。*进多出少，则$theta$增。反之，则降。*
    2. $Delta z_i$为正，注意区分坐标系是*向下为正，还是向上为正*。
  ],
)

将公式#[@q_i0]、#[@q_i1]带入#[@d_theta]，得到：

$
  (Delta theta_(i)) / (Delta t) Delta z_i =
  - & [calqRight(i-1, i, n, n+1)] + \
    & [calqRight(i, i+1, n, n+1)]
      - e_i
$

合并同类项，整理成$a Delta theta_(i-1) + b Delta theta_(i) + c Delta theta_(i+1) = r$的形式：

$
  - pdv(q_(i-1), theta_(i-1)) Delta theta_(i-1) +
  [pdv(q_(i), theta_(i)) - pdv(q_(i-1), theta_(i)) - (Delta z_i) / (Delta t) ] Delta theta_(i) +
  pdv(q_(i), theta_(i+1)) Delta theta_(i+1)
  = q_(i-1)^(n) - q_(i)^(n) + e_i
$

因此，
$
  a & = - pdv(q_(i-1), theta_(i-1)) \
  b & = pdv(q_(i), theta_(i)) - pdv(q_(i-1), theta_(i)) - (Delta z_i) / (Delta t) \
  c & = pdv(q_(i), theta_(i+1)) \
  r & = q_(i-1)^(n) - q_(i)^(n) + e_i
$

之后，可利用带状矩阵的方法进行求解。

两类边界条件。

- *i=1时：*

$
  (Delta theta_(i)) / (Delta t) Delta z_i =
  -q_(i-1)^(n+1) + & [calqRight(i, i+1, n, n+1)]
                     - e_i
$
$
  [pdv(q_i, theta_i) - (Delta z_i) / (Delta t)] Delta theta_(i) +
  pdv(q_i, theta_(i+1)) Delta theta_(i+1) = q_(i-1)^(n+1) - q_i^n + e_i
$

$
  a & = 0, \
  b & = pdv(q_i, theta_i) - (Delta z_i) / (Delta t) \
  c & = pdv(q_i, theta_(i+1)) \
  r & = q_(i-1)^(n+1) - q_i^n + e_i
$
#beamer-block(
  [另外一种理解方式：$q_0$是初始边界条件，不随$theta_0$和$theta_1$变化而变化。],
)

- *i=N时：*

*若底部自由排水*，$q_i^(n+1) = -K_i^(n+1) = -K_I^(n) + pdv(K_i, theta_i) Delta theta_i$，$pdv(q_i^(n+1), theta_i) = -pdv(K_i, theta_i)$、$pdv(q_i^(n+1), theta_(i+1)) = 0$，可得，

$
  a &= - pdv(q_(i-1), theta_(i-1)) \
  b &= pdv(q_(i), theta_(i)) - pdv(q_(i-1), theta_(i)) - (Delta z_i) / (Delta t) = -pdv(K_i, theta_i) - pdv(q_(i-1), theta_(i)) - (Delta z_i) / (Delta t)\
  c &= 0 \
  r &= q_(i-1)^(n) - q_(i)^(n) + e_i
$

*若底部零排水*，$q_i^(n+1) = C = 0$

$
  (Delta theta_(i)) / (Delta t) Delta z_i =
  -[calqRight(i-1, i, n, n+1)] +
  q_i^(n+1) - e_i
$

$
  -pdv(q_(i-1), theta_(i-1)) Delta theta_(i-1)
  -[(Delta z_i) / (Delta t) + pdv(q_(i-1), theta_i) ]Delta theta_i = q_(i-1)^n - q_i^(n+1) + e_i
$

$
  a & = -pdv(q_(i-1), theta_(i-1)), \
  b & = -[(Delta z_i) / (Delta t) + pdv(q_(i-1), theta_i)], \
  c & = 0, \
  r & = q_(i-1)^n - q_i^(n+1) + e_i
$

- *i=N+1时（地下水层，Zeng2009方案）：*

在Zeng2009方案中，为了考虑土壤水与地下水的相互作用，增加了第N+1层（地下水层）。该层位于地下水位$z_"wt"$附近。

*几何设置*：
- 第N+1层中心位置：$z_(N+1) = 0.5(z_"wt" + z_N)$（地下水位与第N层中心的中点）
- 第N+1层厚度：$Delta z_"gw" = |z_"wt" - z_N|$（当$"jwt" >= N$时）

*质量守恒方程*：

第N+1层只有来自第N层的入流$q_N$，没有出流（$q_(N+1) = 0$），且不考虑蒸发/蒸腾（$e_(N+1) = 0$）：

$
  (Delta theta_(N+1)) / (Delta t) Delta z_"gw" = -q_N^(n+1) + 0 - 0 = -q_N^(n+1)
$

*线性化*：

$
  q_N^(n+1) = q_N^n + pdv(q_N, theta_N) Delta theta_N + pdv(q_N, theta_(N+1)) Delta theta_(N+1)
$

代入质量守恒方程：

$
  (Delta z_"gw") / (Delta t) Delta theta_(N+1) = -q_N^n - pdv(q_N, theta_N) Delta theta_N - pdv(q_N, theta_(N+1)) Delta theta_(N+1)
$

整理得：

$
  -pdv(q_N, theta_N) Delta theta_N - [pdv(q_N, theta_(N+1)) + (Delta z_"gw") / (Delta t)] Delta theta_(N+1) = q_N^n
$

*三角矩阵系数*：

$
  a_(N+1) & = -pdv(q_N, theta_N) \
  b_(N+1) & = -pdv(q_N, theta_(N+1)) - (Delta z_"gw") / (Delta t) \
  c_(N+1) & = 0 \
  r_(N+1) & = q_N^n
$

// #pagebreak()

=== 1.4 给水度（specific yield）

靠近地下水水位$z_{"wt"}$的土壤层，认为土壤水势达到均衡状态，即：

#mitex(
  `$$
\psi + z = \psi_{sat} + z_{wt}
$$`,
)

$
  s_(y,i) = & theta_("sat",i) - theta(psi_("sat",i) + z_("wt")) \
          = & theta_("sat",i) - theta(psi_(0))
$

通过积分推导（CoLM，公式12.30~12.34），可以得到，第i层土壤$"layer"_i$的产水率$s_(y,i)$（specific yield）。$s_(y,i)$的物理意义为土壤表层（$z=0$）的空气体积百分比。这主要是因为地下水水位的上升，等价于土壤含水量曲线的整体向上平移，即最上层的干旱部分被地下水附近的饱和部分代替。
