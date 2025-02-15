## 基于 van Genuchten 1980 的均衡水势

> 参考Zeng2009 Campbell1974的方案

<https://chatgpt.com/c/67aad86c-8458-8012-85e7-07d3ae77aa16>

**CHATGPT o3-mini-high**

$$
\int_{C - z_{i+1/2}}^{C - z_{i-1/2}} \left[\theta_r + \frac{\theta_s - \theta_r}{(1 + (-\alpha \psi)^n)^m}\right] d\psi, m = 1 - 1/n
$$
**求解过程：**
$$
I_1 = \int_{C-z_{i+1/2}}^{C-z_{i-1/2}} \theta_r\,d\psi
= \theta_r\Bigl[\psi\Bigr]_{C-z_{i+1/2}}^{C-z_{i-1/2}}
= \theta_r\Bigl[(C-z_{i-1/2})-(C-z_{i+1/2})\Bigr]
= \theta_r\,(z_{i+1/2}-z_{i-1/2})\,.
$$

$$
I_2 = (\theta_s-\theta_r)\int_{C-z_{i+1/2}}^{C-z_{i-1/2}} \bigl(1+(-\alpha\psi)^n\bigr)^{\frac{1}{n}-1}\,d\psi\,.
$$

令$f(\psi)=\bigl(1+(-\alpha\psi)^n\bigr)^{\frac{1}{n}-1}$，求$F$，$F' = f(\psi)$

根据高斯超几何函数（Gaussian hypergeometric function）
$$
F(\psi)=\psi\,_2F_1\!\Bigl(\frac{1}{n},\,1;\,1+\frac{1}{n};\,-(\alpha\psi)^n\Bigr),
$$

$$
\frac{d}{d\psi}\Biggl\{\psi\,_2F_1\!\Bigl(\frac{1}{n},\,1;\,1+\frac{1}{n};\,-(\alpha\psi)^n\Bigr)\Biggr\} =\Bigl(1+(\alpha\psi)^n\Bigr)^{\frac{1}{n}-1}\,.
$$

因此，
$$
I_2 = (\theta_s-\theta_r) [F(C-z_{i-1/2})-F(C-z_{i+1/2})]
$$
也既是最终的解析解为：
$$
\boxed{ \begin{aligned} I &= \theta_r\,(z_{i+1/2}-z_{i-1/2}) \\ &\quad + (\theta_s-\theta_r)\Biggl\{ (C-z_{i-1/2})\,_2F_1\!\Bigl(\frac{1}{n},\,1;\,1+\frac{1}{n};\,-\bigl[\alpha(C-z_{i-1/2})\bigr]^n\Bigr) \\ &\qquad\quad - (C-z_{i+1/2})\,_2F_1\!\Bigl(\frac{1}{n},\,1;\,1+\frac{1}{n};\,-\bigl[\alpha(C-z_{i+1/2})\bigr]^n\Bigr) \Biggr\}\,. \end{aligned} }
$$
