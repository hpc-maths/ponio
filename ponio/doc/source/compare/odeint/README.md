# odeint

> solving ODEs in C++.

| Link     | [Documentation](https://www.boost.org/doc/libs/1_87_0/libs/numeric/odeint) |
|----------|----------------------------------------------------------------------------|
| Language | C++                                                                        |
| License  | Boost Software License                                                     |

## Lorenz equations

We would like to solve the Lorenz equations, a classical chaotic system given by:

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - x z \\
    \dot{z} &= x y - \beta z
  \end{cases}
$$

with parameter $\sigma=10$, $\rho = 28$ and $\beta = \frac{8}{3}$, and the initial state $(x_0, y_0, z_0) = (1,1,1)$. We solve it with a classical Runge-Kutta method of order 4, given by `boost::numeric::odeint::runge_kutta4<state_t>` method.

To solve a problem with odeint we first write our equation as a ODE of the form:

$$
  \dot{u} = f(t, u)
$$

and the user provide the function $f$ as (a lambda function for the example):

```cpp
  auto f = [](state_t const& u, state_t & du, double t);
```

where `u` the current state of the function, `t` the current time and output `du` $f(t,u)$ by reference.

```{literalinclude} lorenz.cpp
  :lines: 15-20
  :language: cpp
  :linenos:
  :lineno-start: 15
```

After defined a method with `boost::numeric::odeint::runge_kutta4<state_t>()`, we can solve the problem between initial time and final time and give an observer which be call after each succeed time iteration

```{literalinclude} lorenz.cpp
  :lines: 37-38
  :language: cpp
  :linenos:
  :lineno-start: 37
```

the `vec_observer` is a lambda function which store all iteration in a `std::vector`.

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

We choose a first order up-wind scheme to estimate the $x$ derivative and a forward Euler method for the time discretization given by `boost::numeric::odeint::euler<state_t>` method.

We define the up-wind scheme as:

```{literalinclude} transport.cpp
  :lines: 48-58
  :language: cpp
  :linenos:
  :lineno-start: 48
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
  :lines: 13-27
  :language: cpp
  :linenos:
  :lineno-start: 13
```

We solve this example with given method `boost::numeric::odeint::runge_kutta_dopri5<state_t>` which is the method RK5(4) 7M in [[DP80](https://doi.org/10.1016/0771-050X(80)90013-3)] (mainly call *DOPRI5*), and need to embedded it into `boost::numeric::odeint::make_controlled` to make an adaptive time step method, and solve the system with a specific function for adaptive time step method:

```{literalinclude} arenstorf.cpp
  :lines: 44-45
  :language: cpp
  :linenos:
  :lineno-start: 44
```

For the complet example, see [`arenstorf.cpp` source file](arenstorf.cpp).
