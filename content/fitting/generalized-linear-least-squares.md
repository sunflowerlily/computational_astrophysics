# General Linear Least Squares

线性最小二乘法中的**线性**指的是参数在拟合函数$Y$中的呈现方式。所以形如：

$$Y(x; \{a_j\}) = \sum_{j=1}^M a_j \varphi_j(x)$$

的式子，即便**基函数** $\{\varphi_j\}$ 是非线性的，但对于 $\{a_j\}$ 而言仍是线性的。


```{admonition} Example
The fitting function:

$$Y(x; \{a_j\}) = a_1 + a_2 x + a_3 x^2$$

is linear in the fit parameters $\{ a_j\}$,  The basis functions in this case are:

\begin{align*}
\varphi_1 &= 1 \\
\varphi_2 &= x \\
\varphi_3 &= x^2
\end{align*}

We can apply the same technique we just did for fitting to a line for this general case.
```


```{admonition} Reference
The discussion in Garcia, _Numerical Methods for Physics_, gives a nice overview, which we loosely follow here.
```

Our $\chi^2$ is:

\begin{align*}
\chi^2(\{a_j\}) &= \sum_{i=1}^N \frac{(Y(x_i; \{a_j\}) - y_i)^2}{\sigma_i^2} \\
    &= \sum_{i=1}^N \frac{1}{\sigma_i^2} \left [
          \left (\sum_{j=1}^M a_j \varphi_j(x_i)\right ) - y_i \right ]^2
\end{align*}

We can differentiate it with respect to one of the parameters, $a_k$:

\begin{align*}
\frac{\partial \chi^2}{\partial a_k}
    &= \frac{\partial}{\partial a_k}
          \sum_{i=1}^N \frac{1}{\sigma_i^2} \left [\left (\sum_{j=1}^M a_j \varphi_j(x_i)\right ) - y_i \right ]^2 \\
    &= \sum_{i=1}^N \frac{1}{\sigma_i^2}
          \frac{\partial}{\partial a_k} \left [\left (\sum_{j=1}^M a_j \varphi_j(x_i)\right ) - y_i \right ]^2 \\
    &= 2 \sum_{i=1}^N \frac{1}{\sigma_i^2} \left [\left (\sum_{j=1}^M a_j \varphi_j(x_i)\right ) - y_i \right ] \varphi_k(x_i) = 0
\end{align*}

or rewriting:

$$\sum_{i=1}^N \sum_{j=1}^M a_j \frac{\varphi_j(x_i) \varphi_k(x_i)}{\sigma_i^2} =
   \sum_{i=1}^N \frac{y_i \varphi_k(x_i)}{\sigma_i^2}$$

We define the $N\times M$ [_design matrix_](https://en.wikipedia.org/wiki/Design_matrix) as

$$A_{ij} = \frac{\varphi_j(x_i)}{\sigma_i}$$

and the source as:

$$b_i = \frac{y_i}{\sigma_i}$$

our system is:

$$\sum_{i=1}^N \sum_{j=1}^M A_{ik} A_{ij} a_j = \sum_{i=1}^N A_{ik} b_i$$

which, by looking at which indices contract, gives us row $k$ of the linear system:

$${\bf A}^\intercal {\bf A} {\bf a} = {\bf A}^\intercal {\bf b}$$

where ${\bf A}^\intercal {\bf A}$ is an $M\times M$ matrix.

The procedure we described above is sometimes called [_ordinary least
squares_](https://en.wikipedia.org/wiki/Ordinary_least_squares).


## Linear fit revisited

For a linear fit,

$$Y(x) = a_1 + a_2 x$$

and our basis functions are: $\phi_1 = 1$ and $\phi_2 = x$.

### Design matrix and source

Our design matrix in this case is:

$${\bf A} = \left ( \begin{array}{cc}
                1/\sigma_1 & x_1 / \sigma_1 \\
                1/\sigma_2 & x_2 / \sigma_2 \\
                \vdots & \vdots \\
                1/\sigma_N & x_N / \sigma_N \\
               \end{array}\right )$$

and the source is:

$${\bf b} = \left (\begin{array}{c} y_1 / \sigma_1 \\
                                    y_2 / \sigma_2 \\
                                    \vdots \\
                                    y_N / \sigma_N \end{array} \right )$$


### Linear system

现在我们可以使用这个设计矩阵和源向量来求解底层的线性系统。

${\bf A}^\intercal {\bf A}$ is:

\begin{align*}
{\bf A}^\intercal{\bf A} &= \left ( \begin{array}{cccc}
                            1/\sigma_1 & 1/\sigma_2 & \cdots & 1/\sigma_N \\
                            x_1/\sigma_1 & x_2/\sigma_2 & \cdots & x_N/\sigma_N \end{array} \right )
                            \left ( \begin{array}{cc}
                1/\sigma_1 & x_1 / \sigma_1 \\
                1/\sigma_2 & x_2 / \sigma_2 \\
                \vdots & \vdots \\
                1/\sigma_N & x_N / \sigma_N \\
               \end{array}\right ) \\
               &= \left ( \begin{array}{cc} \sum_i 1/\sigma_i^2 & \sum_i x_i / \sigma_i^2 \\
                                           \sum_i x_i/\sigma_i^2 & \sum_i x_i^2 / \sigma_i^2 \end{array} \right )
\end{align*}

and ${\bf A}^\intercal {\bf A} {\bf a}$ is:

\begin{align*}
{\bf A}^\intercal {\bf A} {\bf a} &=
   \left ( \begin{array}{cc} \sum_i 1/\sigma_i^2 & \sum_i x_i / \sigma_i^2 \\
                             \sum_i x_i/\sigma_i^2 & \sum_i x_i^2 / \sigma_i^2 \end{array} \right )
   \left ( \begin{array}{c} a_1 \\ a_2 \end{array} \right ) \\
   &= \left ( \begin{array}{c} a_1 \sum_i 1/\sigma_i^2 + a_2 \sum_i x_i/\sigma_i^2 \\
                               a_1 \sum_i x_i/\sigma_i^2 + a_2 \sum_i x_i^2 /\sigma_i^2 \end{array} \right )
                               \end{align*}


${\bf A}^\intercal {\bf b}$ is:

\begin{align*}
{\bf A}^\intercal {\bf b} &= \left ( \begin{array}{cccc}
                            1/\sigma_1 & 1/\sigma_2 & \cdots & 1/\sigma_N \\
                            x_1/\sigma_1 & x_2/\sigma_2 & \cdots & x_N/\sigma_N \end{array} \right )
           \left ( \begin{array}{c}
                y_1 / \sigma_1 \\
                y_2 / \sigma_2 \\
                \vdots \\
                y_N / \sigma_N \\
               \end{array}\right ) \\
               &= \left ( \begin{array}{c} \sum_i y_i / \sigma_i^2 \\
                                           \sum_i x_i y_i / \sigma_i^2 \end{array} \right )
\end{align*}

将这些结合起来，我们得到一个包含两个未知数 ($a_1$, $a_2$) 的二元方程组：

\begin{align*}
a_1 \sum_i \frac{1}{\sigma_i^2} + a_2 \sum_i \frac{x_i}{\sigma_i^2} &= \sum_i \frac{y_i}{\sigma_i^2} \\
a_1 \sum_i \frac{x_i}{\sigma_i^2} + a_2 \sum_i \frac{x_i^2}{\sigma_i^2} &= \sum_i \frac{x_i y_i}{\sigma_i^2}
\end{align*}

这正是我们之前看到的方程组。
