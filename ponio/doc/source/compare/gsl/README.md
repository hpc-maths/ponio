# GSL - GNU Scientific Library

> The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free software under the GNU General Public License.

| Link           | [Web site](https://www.gnu.org/software/gsl/)              |
|----------------|------------------------------------------------------------|
| Language       | C                                                          |
| License        | GPL                                                        |
| Documentation  | [documentation](https://www.gnu.org/software/gsl/doc/html) |

## Lorenz equations

We would like to solve the Lorenz equations, a classical chaotic system given by:

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - x z \\
    \dot{z} &= x y - \beta z
  \end{cases}
$$

with parameter $\sigma=10$, $\rho = 28$ and $\beta = \frac{8}{3}$, and the initial state $(x_0, y_0, z_0) = (1,1,1)$. We solve it with a classical Runge-Kutta method of order 4, given by `gsl_odeiv2_step_rk4` method.

To solve a problem with GSL we first write our equation as a ODE of the form:

$$
  \dot{u} = f(t, u)
$$

and the user provide the function $f$ as:

```c
  int f(double t, const double u[], double du[], void* params);
```

where `u` the current state of the function, `t` the current time and output `du` $f(t,u)$. To use GSL you should convert your state type in C-style array.

```{literalinclude} lorenz.c
  :lines: 7-18
  :language: c
  :linenos:
  :lineno-start: 7
```

Now we should define the system and the driver for the chosen integrator

```{literalinclude} lorenz.c
  :lines: 29-30
  :language: c
  :linenos:
  :lineno-start: 29
```

now we can call the integrator by

```{literalinclude} lorenz.c
  :lines: 53-53
  :language: c
  :linenos:
  :lineno-start: 53
```

The complet time loop is

```{literalinclude} lorenz.c
  :lines: 45-61
  :language: c
  :linenos:
  :lineno-start: 45
```

where `sol` is an array where we save all states.

For the complet example, see [`lorenz.c` source file](lorenz.c).

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

To define the explicit Euler method, we define a ``euler_state_t`` data structure to store current state and intermediate step:

```{literalinclude} transport.c
  :lines: 14-19
  :language: c
```

and add a new ``gsl_odeiv2_step_type`` with following lines:

```{literalinclude} transport.c
  :lines: 157-167
  :language: c
```

where we should define how to allocate a ``euler_state_t`` (``euler_alloc``), increment it with ``euler_apply``, reset data with ``euler_reset``, a function to return the order ``euler_order`` and a deallocator ``euler_free``.

We define the up-wind scheme as:

```{literalinclude} transport.c
  :lines: 169-186
  :language: c
  :linenos:
  :lineno-start: 169
```

The time loop is the same as for Lorenz equation.

For the complet example, see [`transport.c` source file](transport.c).

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

```{literalinclude} arenstorf.c
  :lines: 8-26
  :language: c
  :linenos:
  :lineno-start: 8
```

We solve this example with given method `gsl_odeiv2_step_rk8pd` which is the method RK8(7) 13M in [[PD81](https://doi.org/10.1016/0771-050X(81)90010-3)] (mainly call *DOPRI8*). To do this we define our system to GSL and also a stepper, and some objects to control adaptive time step method:

```{literalinclude} arenstorf.c
  :lines: 37-42
  :language: c
  :linenos:
  :lineno-start: 37
```

then the time loop becomes

```{literalinclude} arenstorf.c
  :lines: 57-75
  :language: cpp
  :linenos:
  :lineno-start: 57
```

we should call `gsl_odeiv2_evolve_apply` to only make a step and not call the method from initial time to final time.

For the complet example, see [`arenstorf.c` source file](arenstorf.c).
