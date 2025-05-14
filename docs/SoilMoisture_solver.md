## 1. 隐式积分

$$
Q_i^{n+1} = -\frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} \left( \psi_i^{n+1} - \psi_{i+1}^{n+1} \right) - K_{i+1/2}^{n+1} \\
Q_{i-1}^{n+1} = -\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} \left( \psi_{i-1}^{n+1} - \psi_i^{n+1} \right) - K_{i-1/2}^{n+1} \\
\frac{\theta_i^{n+1} - \theta_i^n}{\Delta t} = - \frac{Q_{i-1}^{n+1} - Q_i^{n+1}}{\Delta z_i} \\
\theta_i^{n+1} = C_i^{n+1} \psi_i^{n+1}, \quad \theta_i^n = C_i^{n+1} \psi_i^n
$$

联立可得：
$$
\frac{C_i^{n+1} \psi_i^{n+1} - C_i^n \psi_i^n}{\Delta t} = 
- \frac{-\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} \left( \psi_{i-1}^{n+1} 
- \psi_i^{n+1} \right) - K_{i-1/2}^{n+1} + \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} \left( \psi_i^{n+1} - \psi_{i+1}^{n+1} \right) + K_{i+1/2}^{n+1}}{\Delta z_i}
$$

$$
(-C_i^{n+1} \psi_i^{n+1} + C_i^n \psi_i^n)   \frac{\Delta z_i}{\Delta t} = -\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} \left( \psi_{i-1}^{n+1} - \psi_i^{n+1} \right) - K_{i-1/2}^{n+1} + \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} \left( \psi_i^{n+1} - \psi_{i+1}^{n+1} \right) + K_{i+1/2}^{n+1}
$$

联立公式，写成$a \psi_{i-1}^{n+1} + b \psi_i^{n+1} + c \psi_{i+1}^{n+1} = d$的形式
$$
-\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} \psi_{i-1}^{n+1} + 
(
  \frac{C_i^{n+1} \Delta z_i}{\Delta t} + 
  \frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} + 
  \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}}
) \psi_{i}^{n+1} 
- \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} \psi_{i+1}^{n+1} = \\

\frac{C_i^n \Delta z_i} {\Delta t} \psi_i^n + K_{i-1/2}^{n+1} - K_{i+1/2}^{n+1}
$$

## 2. Crank-Nicolson

n+1时刻：

$$
Q_i^{n+1} = -\frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} \left( \psi_i^{n+1} - \psi_{i+1}^{n+1} \right) - K_{i+1/2}^{n+1} \\
Q_{i-1}^{n+1} = -\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} \left( \psi_{i-1}^{n+1} - \psi_i^{n+1} \right) - K_{i-1/2}^{n+1} \\
$$

n时刻：
$$
Q_i^{n} = -\frac{K_{i+1/2}^{n}}{\Delta z_{i+1/2}} \left( \psi_i^{n} - \psi_{i+1}^{n} \right) - K_{i+1/2}^{n} \\
Q_{i-1}^{n} = -\frac{K_{i-1/2}^{n}}{\Delta z_{i-1/2}} \left( \psi_{i-1}^{n} - \psi_i^{n} \right) - K_{i-1/2}^{n} \\
$$

$$
\theta_i^{n+1} = C_i^{n+1} \psi_i^{n+1}, 
\theta_i^n = C_i^{n+1} \psi_i^n
$$

联立公式，

$$
Q_{i-1}^{n+1} - Q_i^{n+1} = \left( -\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} (\psi_{i-1}^{n+1} - \psi_i^{n+1}) - K_{i-1/2}^{n+1} \right) - \left( - \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} (\psi_i^{n+1} - \psi_{i+1}^{n+1}) - K_{i+1/2}^{n+1} \right) \\ 
 = -\frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} (\psi_{i-1}^{n+1} - \psi_i^{n+1}) + \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} (\psi_i^{n+1} - \psi_{i+1}^{n+1}) - (K_{i-1/2}^{n+1} - K_{i+1/2}^{n+1})
$$

$$
Q_{i-1}^{n} - Q_i^{n} = - \frac{K_{i-1/2}^{n}}{\Delta z_{i-1/2}} (\psi_{i-1}^{n} - \psi_i^{n}) + \frac{K_{i+1/2}^{n}}{\Delta z_{i+1/2}} (\psi_i^{n} - \psi_{i+1}^{n}) - (K_{i-1/2}^{n} - K_{i+1/2}^{n})
$$

Crank-Nicolson方式，进行积分：
$$
\frac{\theta_i^{n+1} - \theta_i^n}{\Delta t} = 
  - \frac{  Q_{i-1}^{n+1} - Q_i^{n+1}}{2\Delta z_i} 
  - \frac{  Q_{i-1}^{n} - Q_i^{n}}{2\Delta z_i} 
$$

$$
\frac{-\theta_i^{n+1} + \theta_i^n}{\Delta t} = 
  \frac{  Q_{i-1}^{n+1} - Q_i^{n+1}}{2\Delta z_i} 
  +\frac{  Q_{i-1}^{n} - Q_i^{n}}{2\Delta z_i} 
$$

$$
\frac{-C_i^{n+1} \psi_i^{n+1} + C_i^n \psi_i^n}{\Delta t} = 
\frac{1}{2 \Delta z_i} \left[ \left( - \frac{K_{i-1/2}^{n+1}}{\Delta z_{i-1/2}} (\psi_{i-1}^{n+1} - \psi_i^{n+1}) + \frac{K_{i+1/2}^{n+1}}{\Delta z_{i+1/2}} (\psi_i^{n+1} - \psi_{i+1}^{n+1}) - (K_{i-1/2}^{n+1} - K_{i+1/2}^{n+1}) \right) \right] +  \\

\frac{1}{2 \Delta z_i} \left[ \left( - \frac{K_{i-1/2}^{n}}{\Delta z_{i-1/2}} (\psi_{i-1}^{n} - \psi_i^{n}) + \frac{K_{i+1/2}^{n}}{\Delta z_{i+1/2}} (\psi_i^{n} - \psi_{i+1}^{n}) - (K_{i-1/2}^{n} - K_{i+1/2}^{n}) \right) \right]
$$


写成$a \psi_{i-1}^{n+1} + b \psi_i^{n+1} + c \psi_{i+1}^{n+1} = d$的形式，

$$
 (- \frac{K_{i-1/2}^{n+1}}{2 \Delta z_{i-1/2}} )   \psi_{i-1}^{n+1} + 
 ( \frac{K_{i-1/2}^{n+1}}{2 \Delta z_{i-1/2}} + \frac{K_{i+1/2}^{n+1}}{2 \Delta z_{i+1/2}} + \frac{C_i^{n+1} \Delta z_i}{\Delta t} ) \psi_{i}^{n+1} +
  (- \frac{K_{i+1/2}^{n+1}}{2 \Delta z_{i+1/2}} )   \psi_{i+1}^{n+1} = \\ 

  \frac{C_i^{n} \Delta z_i}{\Delta t} \psi_i^{n} + 
    \frac{K_{i-1/2}^{n}}{2 \Delta z_{i-1/2}} (\psi_{i-1}^{n} - \psi_i^{n}) - 
    \frac{K_{i+1/2}^{n}}{2 \Delta z_{i+1/2}} (\psi_i^{n} - \psi_{i+1}^{n}) +
  (K_{i-1/2}^{n} - K_{i+1/2}^{n}) / 2 + 
  (K_{i-1/2}^{n+1} - K_{i+1/2}^{n+1}) / 2
$$

## 3. Picard iteration

采用第m次的结果，更新K, C，然后重新求解$\psi^{n+1}$
$$
- \frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} \psi^{n+1,m+1}_{i-1} + 
\left( 
	\frac{C_i^{n+1,m} \Delta z_i}{\Delta t} + 
	\frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} + 
	\frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} 
\right) \psi^{n+1,m+1}_{i} 
	- \frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} \psi^{n+1,m+1}_{i+1} = \\
 C_i^{n+1,m} \frac{\Delta z_i}{\Delta t} \psi^n_i + 
 K^{n+1,m}_{i-1/2}
 - K^{n+1,m}_{i+1/2} 
$$

令$\delta^{m+1} = \psi^{n+1,m+1} - \psi^{n+1,m}$，

$$
- \frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} (
  \psi^{n+1, m}_{i-1} + \delta^{m+1}_{i-1}
) + 
\left( 
	\frac{C_i^{n+1,m} \Delta z_i}{\Delta t} + 
	\frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} + 
	\frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} 
\right) (
  \psi^{n+1, m}_i + \delta^{m+1}_{i}
) 
	- \frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} (
    \psi^{n+1, m}_{i+1} + \delta^{m+1}_{i+1}
  ) = \\
 C_i^{n+1,m} \frac{\Delta z_i}{\Delta t} \psi^n_i + 
 K^{n+1,m}_{i-1/2}
 - K^{n+1,m}_{i+1/2} 
$$

展开之后，可得：

$$
- \frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} \delta^{m+1}_{i-1} + 
\left( 
	\frac{C_i^{n+1,m} \Delta z_i}{\Delta t} + 
	\frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} + 
	\frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} 
\right) \delta^{m+1}_{i} 
	- \frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} \delta^{m+1}_{i+1} = \\

  \frac{K^{n+1,m}_{i-1/2}}{\Delta z_{i-1/2}} (\psi^{n+1, m}_{i-1} - \psi^{n+1, m}_{i}) - 
  \frac{K^{n+1,m}_{i+1/2}}{\Delta z_{i+1/2}} (\psi^{n+1, m}_{i} - \psi^{n+1, m}_{i+1})  \\

 C_i^{n+1,m} \frac{\Delta z_i}{\Delta t} (\psi^n_i - \psi^{n+1,m}_i) + 
 K^{n+1,m}_{i-1/2}
 - K^{n+1,m}_{i+1/2}
$$
