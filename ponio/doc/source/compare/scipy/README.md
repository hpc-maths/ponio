# SciPy

> Fundamental algorithms for scientific computing in Python.

| Link           | [Web site](https://scipy.org/)                                                   |
|----------------|----------------------------------------------------------------------------------|
| Language       | Python                                                                           |
| License        | BSD-3-Clause license                                                             |
| Documentation  | [documentation](https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide) |
| Git repository | [@scipy](https://github.com/scipy/scipy)                                         |

## Lorenz equations
We would like to solve the Lorenz equations, a classical chaotic system given by:

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - x z \\
    \dot{z} &= x y - \beta z
  \end{cases}
$$

with parameter $\sigma=10$, $\rho = 28$ and $\beta = \frac{8}{3}$, and the initial state $(x_0, y_0, z_0) = (1,1,1)$. We solve it with a classical Runge-Kutta method of order 4.

To solve Lorenz equations with a classical 4th order Runge-Kutta RK(4,4) we need to define it:

```{literalinclude} lorenz.py
  :lines: 6-17
```

and now call it with `scipy.integrate.solve_ivp` function:

```python
sol = solve_ivp(lorenz_system, t_span, y0, method=RK44, first_step=dt)
```

To solve a problem with SciPy we first write our equation as a ODE of the form:

$$
  \dot{u} = f(t, u)
$$

and the user provide the function $f$ as:

```py
  def f(t, u):
    # ...
    return du
```

where `u` the current state of the function, `t` the current time and output `du` $f(t,u)$.

```{literalinclude} lorenz.py
  :lines: 27-36
  :language: py
  :linenos:
  :lineno-start: 27
```

After chose a method, we can solve the problem between initial time and final time

```{literalinclude} lorenz.py
  :lines: 43-44
  :language: py
  :linenos:
  :lineno-start: 43
```

this function returns an object with solution at each time (and also dense output properties).

For the complet example, see [`lorenz.py` source file](lorenz.py).

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

We choose a first order up-wind scheme to estimate the $x$ derivative and a forward Euler method for the time discretization.

SciPy provides mainly adaptive time step methods, so there is no classical 4th order Runge-Kutta RK(4,4) nor explicit Euler (or forward Euler or RK(1, 1)), we need to define it:

```{literalinclude} transport.py
  :lines: 6-13
```

and now call it with `scipy.integrate.solve_ivp` function:

```python
sol = solve_ivp(upwind, t_span, y0, method=Euler, first_step=dt)
```

We define the up-wind scheme as:

```{literalinclude} transport.py
  :lines: 38-48
  :language: py
  :linenos:
  :lineno-start: 38
```

The time loop is the same as for Lorenz equation.

For the complet example, see [`transport.py` source file](transport.py).

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

```{literalinclude} arenstorf.py
  :lines: 6-24
  :language: py
  :linenos:
  :lineno-start: 6
```

We solve this example with given method `RK45` which is the method RK5(4) 7M in [[DP80](https://doi.org/10.1016/0771-050X(80)90013-3)] (mainly call *DOPRI5*) and `DOP853` which is the method RK8(7) 13M in [[PD81](https://doi.org/10.1016/0771-050X(81)90010-3)] (mainly call *DOPRI8*).

The time loop is the same as for Lorenz equation, for `RK45` method

```{literalinclude} arenstorf.py
  :lines: 32-33
  :language: py
  :linenos:
  :lineno-start: 32
```

and for `DOP853` method

```{literalinclude} arenstorf.py
  :lines: 42-43
  :language: py
  :linenos:
  :lineno-start: 42
```

For the complet example, see [`arenstorf.py` source file](arenstorf.py).
