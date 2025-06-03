# ponio

> The library ponio is a collection of time integrators for solving differential equations written in C++

| Link          | [Github repository](https://github.com/hpc-maths/ponio)   |
|---------------|-----------------------------------------------------------|
| Language      | C++                                                       |
| License       | BSD-3-Clause                                              |
| Documentation | [documentation](https://ponio.readthedocs.io/en/latest/)  |
| Analysis      | [analysis of methods](https://hpc-maths.github.io/ponio/) |

## Lorenz equations

We would like to solve the Lorenz equations, a classical chaotic system given by:

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - x z \\
    \dot{z} &= x y - \beta z
  \end{cases}
$$

with parameter $\sigma=10$, $\rho = 28$ and $\beta = \frac{8}{3}$, and the initial state $(x_0, y_0, z_0) = (1,1,1)$. We solve it with a classical Runge-Kutta method of order 4, given by `ponio::runge_kutta::rk_44` method.

To solve a problem with ponio we first write our equation as a ODE of the form:

$$
  \dot{u} = f(t, u)
$$

and the user provide the function $f$ as (a lambda function for the example):

```cpp
  auto f = [](double t, state_t const& u) -> state_t;
```

where `u` the current state of the function, `t` the current time and output is the result of $f(t,u)$.

```{literalinclude} lorenz.cpp
  :lines: 19-22
  :language: cpp
  :linenos:
  :lineno-start: 19
```

After defined a method with `ponio::runge_kutta::rk_44()`, we can solve the problem between initial time and final time and give an observer which be call after each succeed time iteration

```{literalinclude} lorenz.cpp
  :lines: 30-30
  :language: cpp
  :linenos:
  :lineno-start: 30
```

the `ponio::observer::file_observer` is an observer that save all states in a file.

For the complet example, see [`lorenz.cpp` source file](lorenz.cpp).

## Transport equation

In this example we would like to solve the following PDE:

$$
  \partial_t u + a \partial_x u = 0
$$

with $t>0$, on the torus $x\in[0, 1)$, a velocity $a=1$ and the initial condition given by a hat function:

$$
  u(0, x) = \begin{cases}
      x - 0.25  & \text{if } x\in[0.25, 0.5[ \\
      -x + 0.75 & \text{if } x\in[0.5, 0.75[ \\
      0         & \text{else}
  \end{cases}
$$

We choose a first order up-wind scheme to estimate the $x$ derivative and a forward Euler method for the time discretization given by `ponio::runge_kutta::euler` method.

We define the up-wind scheme as:

```{literalinclude} transport.cpp
  :lines: 50-64
  :language: cpp
  :linenos:
  :lineno-start: 50
```

The time loop is the same as for Lorenz equation.

For the complet example, see [`transport.cpp` source file](transport.cpp).

## Arenstorf orbit

The Arenstorf orbit problem is a classical problem to test adaptive time step methods:

$$
  \begin{cases}
    \ddot{x} &= x + 2\dot{y} - \frac{1-\mu}{r_1^3}(x+\mu) - \frac{\mu}{r_2^3}(x-1+\mu) \\
    \ddot{y} &= y - 2\dot{x} - \frac{}{1-\mu}{r_1^3}y - \frac{\mu}{r_2^3}y
  \end{cases}
$$

with initial condition $(x,\dot{x},y,\dot{y})=(0.994, 0, 0, -2.001585106)$, $r_1$ and $r_2$ given by

$$
  r_1 = \sqrt{(x+\mu)^2 + y^2},\quad r_2 = \sqrt{(x-1+\mu)^2 + y^2}
$$

and with parameter $\mu = 0.012277471$.

First of all, we need to rewrite this problem into a first order derivative equation in time

$$
  \begin{pmatrix}
    y_1 \\
    y_2 \\
    y_3 \\
    y_4
  \end{pmatrix}
  =
  \begin{pmatrix}
    x \\
    y \\
    \dot{x} \\
    \dot{y}
  \end{pmatrix},
  \qquad
  \begin{cases}
    \dot{y}_1 = y_3 \\
    \dot{y}_2 = y_4 \\
    \dot{y}_3 = y_1 + 2y_4 - \frac{1-\mu}{r_1^3}(y_1 + \mu) - \frac{\mu}{r_2^3}(y_1-1+\mu) \\
    \dot{y}_4 = y_2 - 2y_3 - \frac{1-\mu}{r_1^3}y_2 - \frac{\mu}{r_2^3}y_2 \\
  \end{cases}
$$

We define this system as:

```{literalinclude} arenstorf.cpp
  :lines: 17-33
  :language: cpp
  :linenos:
  :lineno-start: 17
```

We solve this example with given method `ponio::runge_kutta::rk54_7m` which is the method RK5(4) 7M in [[DP80](https://doi.org/10.1016/0771-050X(80)90013-3)] (mainly call *DOPRI5*) and `ponio::runge_kutta::rk87_13m` which is the method RK8(7) 13M in [[PD81](https://doi.org/10.1016/0771-050X(81)90010-3)] (mainly call *DOPRI8*).

The time loop is the same as for Lorenz equation, for `rk54_7m` method

```{literalinclude} arenstorf.cpp
  :lines: 41
  :language: cpp
  :linenos:
  :lineno-start: 41
```

and for `rk87_13m` method

```{literalinclude} arenstorf.cpp
  :lines: 42
  :language: cpp
  :linenos:
  :lineno-start: 42
```

For the complet example, see [`arenstorf.cpp` source file](arenstorf.cpp).
