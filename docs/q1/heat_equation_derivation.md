# 热传导方程的底层推导：从能量守恒到圆柱坐标热方程

这个热传导方程和水分扩散的逻辑几乎一模一样，只是把“质量守恒 + Fick 定律”换成“能量守恒 + Fourier 定律”。

我们最终使用的控制方程是

\[
\rho c_p\frac{\partial T}{\partial t}
=
\frac1r\frac{\partial}{\partial r}
\left(
kr\frac{\partial T}{\partial r}
\right),
\]

然后因为题目中 \(\rho,c_p,k\) 都取常数，才化成

\[
\frac{\partial T}{\partial t}
=
\alpha
\left(
\frac{\partial^2T}{\partial r^2}
+\frac1r\frac{\partial T}{\partial r}
\right),
\qquad
\alpha=\frac{k}{\rho c_p}.
\]

---

## 1. 第一步：能量守恒

考虑一个很小的体积元。

单位体积内储存的显热，可以写成

\[
e=\rho c_p T.
\]

这里暂时假设 \(\rho\) 和 \(c_p\) 不随温度变化。

于是单位体积能量随时间的变化率为

\[
\frac{\partial e}{\partial t}
=
\rho c_p\frac{\partial T}{\partial t}.
\]

如果内部没有热源，那么这个小体积的能量增加，只能来自：

> 流入的热量大于流出的热量。

因此局部能量守恒写成

\[
\boxed{
\rho c_p\frac{\partial T}{\partial t}
=
-\nabla\cdot\mathbf q
}
\]

其中 \(\mathbf q\) 是热流密度，单位为 W/m\(^2\)。

这一式子和水分扩散中的

\[
\frac{\partial C}{\partial t}
=
-\nabla\cdot\mathbf J
\]

完全对应。

---

## 2. 第二步：Fourier 导热定律

只有守恒方程还不够，因为我们还不知道热流 \(\mathbf q\) 是多少。

Fourier 导热定律给出

\[
\boxed{
\mathbf q=-k\nabla T
}
\]

它表示：热量沿温度降低的方向传递。

如果只看径向，则有

\[
q_r''=-k\frac{\partial T}{\partial r}.
\]

其中：

- \(q_r''\)：径向热通量；
- \(k\)：导热系数；
- \(\partial T/\partial r\)：径向温度梯度。

负号表示热量从高温区流向低温区。

---

## 3. 把 Fourier 定律代入能量守恒

从

\[
\rho c_p\frac{\partial T}{\partial t}
=
-\nabla\cdot\mathbf q
\]

出发，代入

\[
\mathbf q=-k\nabla T,
\]

得到

\[
\rho c_p\frac{\partial T}{\partial t}
=
-\nabla\cdot(-k\nabla T).
\]

两个负号抵消：

\[
\boxed{
\rho c_p\frac{\partial T}{\partial t}
=
\nabla\cdot(k\nabla T)
}
\]

这才是更一般的热传导方程。

它与水分扩散方程

\[
\frac{\partial C}{\partial t}
=
\nabla\cdot(D\nabla C)
\]

具有完全相同的数学结构。

两者可以对应为

\[
C
\leftrightarrow T,
\qquad
D\leftrightarrow k,
\qquad
\mathbf J\leftrightarrow\mathbf q.
\]

区别在于热传导方程左侧多了

\[
\rho c_p,
\]

因为温度 \(T\) 本身不是能量密度。单位温度变化对应多少储存能量，要乘上体积热容 \(\rho c_p\)。

---

## 4. 从三维形式变成圆柱径向形式

一般三维热传导方程为

\[
\rho c_p\frac{\partial T}{\partial t}
=
\nabla\cdot(k\nabla T).
\]

对于圆柱形药材，我们假设：

1. 关于中心轴轴对称；
2. 忽略 \(\theta\) 方向变化；
3. 忽略轴向 \(z\) 方向变化；
4. 温度只依赖于 \(r\) 和 \(t\)。

因此

\[
T=T(r,t).
\]

圆柱坐标下，如果只有径向热流，散度为

\[
\nabla\cdot\mathbf q
=
\frac1r\frac{\partial}{\partial r}(rq_r).
\]

于是能量守恒变成

\[
\rho c_p\frac{\partial T}{\partial t}
=
-\frac1r\frac{\partial}{\partial r}(rq_r).
\]

再代入 Fourier 定律

\[
q_r=-k\frac{\partial T}{\partial r},
\]

得到

\[
\rho c_p\frac{\partial T}{\partial t}
=
-\frac1r
\frac{\partial}{\partial r}
\left[
r\left(-k\frac{\partial T}{\partial r}\right)
\right].
\]

两个负号抵消：

\[
\boxed{
\rho c_p\frac{\partial T}{\partial t}
=
\frac1r\frac{\partial}{\partial r}
\left(
kr\frac{\partial T}{\partial r}
\right)
}
\]

这就是轴对称圆柱中的瞬态热传导方程。

---

## 5. 为什么最后能变成经典热方程

题目中把

\[
k,\rho,c_p
\]

都视为常数。

因此 \(k\) 可以从空间求导中提出：

\[
\rho c_p T_t
=
k\frac1r\frac{\partial}{\partial r}
\left(
rT_r
\right).
\]

而

\[
\frac1r\frac{\partial}{\partial r}(rT_r)
=
T_{rr}+\frac1rT_r.
\]

因此

\[
\rho c_p T_t
=
k
\left(
T_{rr}+\frac1rT_r
\right).
\]

两边除以 \(\rho c_p\)：

\[
T_t
=
\frac{k}{\rho c_p}
\left(
T_{rr}+\frac1rT_r
\right).
\]

定义热扩散率

\[
\boxed{
\alpha=\frac{k}{\rho c_p}
}
\]

于是得到

\[
\boxed{
\frac{\partial T}{\partial t}
=
\alpha
\left(
\frac{\partial^2T}{\partial r^2}
+\frac1r\frac{\partial T}{\partial r}
\right)
}
\]

这就是经典的圆柱坐标热传导方程。

---

## 6. 和水分扩散的对应关系

热传导：

\[
\boxed{
\text{能量守恒}
}
\]

\[
\rho c_p T_t=-\nabla\cdot\mathbf q
\]

加上

\[
\boxed{
\text{Fourier 定律}
}
\]

\[
\mathbf q=-k\nabla T
\]

得到

\[
\boxed{
\rho c_pT_t=\nabla\cdot(k\nabla T)
}
\]

如果 \(k,\rho,c_p\) 都是常数，

\[
\boxed{
T_t=\alpha\nabla^2T
}
\]

其中

\[
\alpha=\frac{k}{\rho c_p}.
\]

水分扩散则是

\[
\frac{\partial C}{\partial t}
=
-\nabla\cdot\mathbf J,
\]

\[
\mathbf J=-D\nabla C,
\]

所以

\[
\boxed{
\frac{\partial C}{\partial t}
=
\nabla\cdot(D\nabla C)
}
\]

因此两者本质上都是：

\[
\boxed{
\text{守恒定律}
+
\text{梯度驱动的通量定律}
}
\]

---

## 7. 如果 \(k\) 不是常数怎么办？

如果以后导热系数变成

\[
k=k(T),
\]

那么就和水分扩散中

\[
D=D(C)
\]

完全一样。

这时不能简单写成

\[
T_t
=
\frac{k(T)}{\rho c_p}\nabla^2T.
\]

正确形式应该保留为

\[
\boxed{
\rho c_pT_t
=
\nabla\cdot[k(T)\nabla T]
}
\]

因为 \(k(T)\) 会通过

\[
T=T(r,t)
\]

间接随空间和时间变化。

也就是说：

- 常 \(k\)：可以化成经典的 \(\alpha\nabla^2T\)；
- 变 \(k\)：必须保留在散度算子内部。

这和变扩散系数 \(D(C)\) 的处理逻辑完全一致。

---

## 总结

热传导方程的底层逻辑可以写成

\[
\boxed{
\text{能量守恒}
}
\]

\[
\Downarrow
\]

\[
\rho c_pT_t=-\nabla\cdot\mathbf q
\]

\[
\Downarrow
\]

\[
\boxed{
\text{Fourier 定律：}
\mathbf q=-k\nabla T
}
\]

\[
\Downarrow
\]

\[
\boxed{
\rho c_pT_t=\nabla\cdot(k\nabla T)
}
\]

\[
\Downarrow
\qquad
(k,\rho,c_p=\text{常数})
\]

\[
\boxed{
T_t=\alpha\nabla^2T
}
\]

其中

\[
\boxed{
\alpha=\frac{k}{\rho c_p}
}
\]

所以经典热方程本身不是最底层的出发点，而是“能量守恒 + Fourier 导热定律”在常物性条件下得到的特殊形式。
