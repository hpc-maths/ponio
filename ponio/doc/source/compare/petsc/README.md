# PETSc - Portable, Extensible Toolkit for Scientific Computation

> PETSc, pronounced PET-see (the S is silent), is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations.

| Link           | [Web site](https://petsc.org/release/)   |
|----------------|------------------------------------------|
| Language       | C                                        |
| License        |  BSD 2-clause license                    |
| Repository     | [GitLab](https://gitlab.com/petsc/petsc) |

## Lorenz equations

We would like to solve the Lorenz equations, a classical chaotic system given by:

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - x z \\
    \dot{z} &= x y - \beta z
  \end{cases}
$$

with parameter $\sigma=10$, $\rho = 28$ and $\beta = \frac{8}{3}$, and the initial state $(x_0, y_0, z_0) = (1,1,1)$. We solve it with a classical Runge-Kutta method of order 4, given by `TSRK4` method.

To solve a problem with PETSc we first write our equation as a ODE of the form:

$$
  \dot{u} = f(t, u)
$$

and the user provide the function $f$ as:

```c
  PetscErrorCode f(TS ts, double t, Vec u, Vec du, void* params);
```

where `u` the current state of the function, `t` the current time and output `du` $f(t,u)$, and `ts` the time stepper.

```{literalinclude} lorenz.c
  :lines: 3-23
  :language: c
  :linenos:
  :lineno-start: 3
```

Next we need to define a time stepper and set all parameters, and give it the $f$ function

```{literalinclude} lorenz.c
  :lines: 43-52
  :language: c
  :linenos:
  :lineno-start: 43
```

now we can call the time stepper by:

```c
  TSStep(ts);
```

and get current state and time with `TSGetSolution` `TSGetTime` functions. The complet time loop is

```{literalinclude} lorenz.c
  :lines: 66-82
  :language: c
  :linenos:
  :lineno-start: 66
```

where `sol` is an array to store all states.

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

We choose a first order up-wind scheme to estimate the $x$ derivative and a forward Euler method for the time discretization given by `TSRK1FE` method.

We define the up-wind scheme as:

```{literalinclude} transport.c
  :lines: 3-29
  :language: c
  :linenos:
  :lineno-start: 3
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
  :lines: 3-30
  :language: c
  :linenos:
  :lineno-start: 3
```

We solve this example with given method `TSRK5DP` which is the method RK8(7) 13M in [[PD81](https://doi.org/10.1016/0771-050X(81)90010-3)] (mainly call *DOPRI8*). For this we need to set the integrator and settings for adaptive time step method

```{literalinclude} arenstorf.c
  :lines: 50-60
  :language: c
  :linenos:
  :lineno-start: 50
```

then the time loop becomes

```{literalinclude} arenstorf.c
  :lines: 74-93
  :language: c
  :linenos:
  :lineno-start: 74
```

For the complet example, see [`arenstorf.c` source file](arenstorf.c).
