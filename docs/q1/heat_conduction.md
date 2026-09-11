# A 题问题 1 热传导模型说明

## 1. 文档目的与建模口径

本文档给出问题 1 前 1800 s 温度场的统一数学模型，作为代码实现、数值验证和论文撰写依据。

问题 1 只考虑药材内部导热和表面对流换热，不考虑蒸发潜热。温度模型与水分模型独立求解，不设置潜热参数、潜热通量或含潜热对照分支。

采用这一简化的原因是题面没有给出真实蒸发质量通量、干物质体积密度、气固界面平衡关系和汽化潜热。水分方程中的 Robin 通量以含水率为状态变量，单位为 m/s，不能直接作为能量方程中的 kg/(m²·s) 质量通量。

## 2. 几何与状态变量

药材近似为固定尺寸圆柱体：

\[
L=0.25\ \mathrm m,
\qquad
R=0.02\ \mathrm m.
\]

建立圆柱坐标系，中心轴为 \(r=0\)，侧表面为 \(r=R\)。忽略周向和轴向变化，只求解

\[
T=T(r,t),
\qquad
0\le r\le R,
\qquad
0\le t\le1800\ \mathrm s.
\]

## 3. 基本假设

1. 药材为均匀、各向同性的有效连续介质。
2. 药材和烘房环境关于圆柱中心轴对称。
3. 温度只沿径向变化，忽略端面换热和轴向导热。
4. 前 1800 s 内圆柱尺寸保持不变，不考虑收缩。
5. 密度、比热容和导热系数取题目给出的常数。
6. 药材内部无体积热源，只考虑傅里叶导热。
7. 侧表面只考虑与烘房空气之间的对流换热。
8. 烘房温度采用附件 1 数据，并在采样点之间作分段线性插值。
9. 不考虑热辐射、反应热、宏观流体运动和蒸发潜热。

## 4. 参数与符号

| 符号 | 含义 | 单位 | 取值或来源 |
| --- | --- | --- | --- |
| \(r\) | 径向坐标 | m | \(0\le r\le0.02\) |
| \(t\) | 时间 | s | \(0\le t\le1800\) |
| \(T(r,t)\) | 药材温度 | °C | 待求 |
| \(T_s(t)\) | 药材表面温度 | °C | \(T(R,t)\) |
| \(T_\infty(t)\) | 烘房空气温度 | °C | 附件 1 |
| \(\rho\) | 药材密度 | kg/m³ | 820 |
| \(c_p\) | 比热容 | J/(kg·K) | 2600 |
| \(k\) | 导热系数 | W/(m·K) | 0.36 |
| \(h\) | 表面对流换热系数 | W/(m²·K) | 25 |
| \(\alpha\) | 热扩散率 | m²/s | \(k/(\rho c_p)\) |

热扩散率为

\[
\alpha
=
\frac{k}{\rho c_p}
\approx1.69\times10^{-7}\ \mathrm{m^2/s}.
\]

## 5. 控制方程

傅里叶定律为

\[
q_r''=-k\frac{\partial T}{\partial r}.
\]

对轴对称薄圆柱壳应用能量守恒，得到

\[
\rho c_p\frac{\partial T}{\partial t}
=
\frac1r\frac{\partial}{\partial r}
\left(
kr\frac{\partial T}{\partial r}
\right).
\]

由于 \(\rho,c_p,k\) 为常数，也可写成

\[
\frac{\partial T}{\partial t}
=
\alpha
\left(
\frac{\partial^2T}{\partial r^2}
+\frac1r\frac{\partial T}{\partial r}
\right).
\]

## 6. 初始条件与边界条件

初始温度均匀：

\[
T(r,0)=28\ ^\circ\mathrm C.
\]

圆柱中心采用对称条件：

\[
\left.\frac{\partial T}{\partial r}\right|_{r=0}=0.
\]

规定径向向外为正方向，侧表面的导热通量与空气对流换热通量平衡：

\[
-k\left.\frac{\partial T}{\partial r}\right|_{r=R}
=
h\left[T_s(t)-T_\infty(t)\right].
\]

等价地，从空气向药材输入热量的方向可写为

\[
k\left.\frac{\partial T}{\partial r}\right|_{r=R}
=
h\left[T_\infty(t)-T_s(t)\right].
\]

## 7. 完整初边值问题

\[
\begin{cases}
\displaystyle
\rho c_p\frac{\partial T}{\partial t}
=
\frac1r\frac{\partial}{\partial r}
\left(kr\frac{\partial T}{\partial r}\right),
&0<r<R,\[1em]
T(r,0)=28,
&0\le r\le R,\[0.5em]
\displaystyle
\left.\frac{\partial T}{\partial r}\right|_{r=0}=0,
&t>0,\[1em]
\displaystyle
-k\left.\frac{\partial T}{\partial r}\right|_{r=R}
=h[T(R,t)-T_\infty(t)],
&t>0.
\end{cases}
\]

## 8. 全局能量守恒

忽略端面后，单位长度圆柱的显热变化应等于侧表面的累计对流输入：

\[
\frac{\mathrm d}{\mathrm dt}
\int_0^R \rho c_pT(r,t)\,2\pi r\,\mathrm dr
=
2\pi Rh[T_\infty(t)-T_s(t)].
\]

离散程序应分别累计储热量和对流输入，并用二者的相对残差检查能量守恒。

## 9. 数值实现

程序采用单元中心有限体积法处理圆柱径向几何，并在 \(r=R\) 处合并半网格导热热阻与表面对流热阻。时间推进采用 Crank–Nicolson 格式，并用两个 Rannacher 半步抑制初始边界变化可能造成的数值振荡。

计算网格和输出网格相互独立：内部使用足够细的径向网格，最终插值到 0–2 cm、间隔 0.1 cm 的官方输出位置，并保存 1–1800 s 的结果。

## 10. 验收标准

- 初始时所有位置均为 \(28\ ^\circ\mathrm C\)。
- 中心热通量为零，表面导热通量与对流换热通量一致。
- 在附件 1 的升温阶段，表面温度通常比中心响应更快。
- 温度场不应出现无依据的负值、过冲或振荡。
- 储热量与累计对流输入满足离散能量守恒。
- 至少比较两组空间网格和两组时间精度设置。
- 恒定环境温度算例应通过圆柱 Bessel 解析解交叉验证。
- 输出工作簿的时间、半径、尺寸、数值类型和保留位数符合官方模板。

## 11. 适用范围

该模型适用于问题 1 的固定尺寸、常物性、短时预热过程。它不描述轴向效应、端面换热、收缩、辐射、材料非均匀性和相变吸热。这些均属于模型简化，应在论文中明确说明，但不在问题 1 代码中增设额外分支。
