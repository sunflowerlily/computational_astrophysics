# 偏微分方程
我们现在将学习偏微分方程。我们将研究 三大类PDEs ：双曲型、椭圆型和抛物型PDEs。每一类都需要不同的求解方法。

```{seealso}
关于天体物理中PDE方法的更深入讨论，可以阅读我的笔记：[_计算天体物理流体力学导论_](https://open-astrophysics-bookshelf.github.io/numerical_exercises/)。
```

## 双曲型PDEs
典型的双曲型PDE是波动方程：

$$\frac{\partial^2 \phi}{\partial t^2} = c^2 \frac{\partial^2 \phi}{\partial x^2}$$

它的解是以速度$c$向两个方向传播的行波：

$$\phi(x, t) = A f_0(x - ct) + B g_0(x + ct)$$

我们将要研究的一个更简单的例子是_线性对流方程_：

$$\frac{\partial a}{\partial t} + u \frac{\partial a}{\partial x} = 0$$

双曲型PDE的定义特征是解的变化具有有限的、真实的传播速度。

在天体物理学中，对流类方程描述了流体动力学。

          
## 椭圆型PDEs

泊松方程是典型的椭圆型PDE：

$$\nabla^2 \phi = f$$

注意这里没有时间变量。关于解的信息以无限快的速度传播，解只依赖于边界条件和源项$f$。

在天体物理学中，泊松方程用于获得引力势，通过[_泊松方程_](https://en.wikipedia.org/wiki/Poisson%27s_equation#Newtonian_gravity)：

$$\nabla^2 \Phi = 4 \pi G \rho$$

## 抛物型PDEs

热传导方程是典型的抛物型PDE：

$$\frac{\partial \phi}{\partial t} = k \frac{\partial^2 \phi}{\partial x^2}$$

与椭圆型PDE类似，关于解的信息是瞬时传递的——不像双曲型那样有有限的传播速度。但这是一个与时间相关的问题。

在恒星中，光子或热的扩散是抛物型PDE在天体物理中出现的主要领域。
        
