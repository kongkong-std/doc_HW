# Numerical integration
Two numerical methods for computing integration, here lists two methods:
1. Trapezoid method
2. Simpson method

Consider 1D integration
$$
\begin{align}
    \mathbf I (f) &=\int _ a ^ b f(x) \ \mathrm d x.
\end{align}
$$


## Trapezoid method
Meshing 1D interval $[a, b]$
$$
\begin{align}
    x _ i &= a + h i, \quad i = 0, 1, \ldots, n,
    \label{eqn:mesh}
\end{align}
$$
where $h = (b - a) / n$, $n$ is the number of elements.

Computing numerical integration $\mathbf I _ n(f)$ with trapezoid method as follows
$$
\begin{align}
    \mathbf I _ n (f) &= \sum _ {i = 0} ^ {n - 1}
    \frac{h}{2} \left[f(x _ i) + f(x _ {i + 1})\right].
\end{align}
$$

## Simpson method
Meshing as previous trapezoid method $\eqref{eqn:mesh}$described. Simpson method as follows
$$
\begin{align}
    \mathbf I _ n (f) &= \sum _ {i = 0} ^ {n - 1}
    h \left[\frac{1}{6} f(x _ i) +
    \frac{4}{6} f\left(\frac{x _ i + x _ {i + 1}}{2}\right) +
    \frac{1}{6} f(x _ {i + 1})\right].
\end{align}
$$

## Task
Integral task 1 is given as follow
$$
\begin{align}
    \mathbf I &= \int _ 0 ^ 1 \frac{\sin x}{x} \ \mathrm d x.
\end{align}
$$
task 2 is described as follow
$$
\begin{align}
    \mathbf I &= \int _ 0 ^ 1 e ^ {x ^ 2} \ \mathrm d x.
\end{align}
$$

## Command line
```linux
./app_exe -n <num_1> -method <string_1> -problem <string_2>
```
where,
1. <num_1> number of elements
2. <string_1> numerical integral method, trapezoid or simpson
3. <string_2> integral problem task_1, task_2