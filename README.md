# ponio

> Presentation Overview of Numerical Integrator for ODE

The library ponio is a collection of time integrators for solving ODE and PDE written in C++. This library aims to be the easiest to use without compromising on performance.

## Installation

### From conda

>  **SOON**

```
  conda install ponio
```

### From source

```
  git clone https://github.com/hpc-maths/ponio.git
  cd ponio
  pixi install
  pixi build
```


## Get started

### Lorenz equations

In this section, we will present how to solve the Lorenz equations, defined:

$$
  \begin{cases}
    \dot{u}_1 &= \sigma(u_2 - u_1) \\
    \dot{u}_2 &= \rho u_1 - u_2 - u_1 u_3 \\
    \dot{u}_3 &= u_1 u_2 - \beta u_3
  \end{cases}
$$

The ponio library solve a problem written of the form:

$$
  \dot{u} = f(t, u)
$$

> The following steps describe how to solve this problem with ponio. It is important to note that these steps are generally the same whatever the equations we want to solve.

First of all, you have to specify the data type that represents a state $u$ of your ODE (generally named `state_t`). This type needs to support arithmetic operations (addition and multiplication by a scalar), that why we use a `std::valarray<double>` in this example.

```cpp
  using state_t = std::valarray<double>;
```

To integrate a differential equation numerically, one also has to define the rhs of the equation $\dot{u} = f(t, u)$. In ponio you supply this function in terms of functor (an object that implements the `()` operator), a function or a lambda function with a certain parameter structure (time and a state). Hence, the straightforward way would be to just define a lambda function, e.g:

```cpp
  const double sigma = 10.;
  const double rho   = 28.;
  const double beta  = 8./3.;

  auto lorenz_rhs = [=]( double /* t */, state_t&& u ) -> state_t
  {
    return {
      sigma * ( u[1] - u[0] ),
      rho * u[0] - u[1] - u[0] * u[2],
      u[0] * u[1] - beta * u[2]
    };
  };
```

Numerical integration works iteratively, that means you start at a state $u(t^n)$ and perform a time-step of length $\Delta t$ to obtain the approximate state $u(t^n+\Delta t) = u(t^{n+1})$. The library ponio calls for you iteratively the method between the initial time $t^0$ the final time $t^N$, which are defined in a time span. You should also define your initial condition $u(t^0)$.

```cpp
  state_t u_ini = {1., 1., 1.};

  double dt = 0.01;
  ponio::time_span<double> t_span = { 0., 20. };
```

> Even if you will use an adaptive time step method to solve the ODE, in ponio you should define an initial time step `dt`.

> The time span can contains intermediate value where the time stepper should pass.

The function `ponio::solve`, which be used to solve the ODE, returns the state at final time, but maybe you want some information on state at each iteration. To do this, ponio library offers some *observers*. For the sake of simplicity we will prestent the file observer which write data in columns, the first one is the current time, the following columns are state (3 columns for the state of Lorenz equation), and the last one is the time step. The ponio library offers literal operators to create a file observer from a string which contains the filename of output file.

```cpp
  using namespace observer;

  auto obs = "sol.txt"_fobs
```

Now you are ready to solve your problem with an explicit Runge-Kutta method, see [algorithm overview](https://ponio.readthedocs.io/en/latest/api/algorithm.html). In this example we will use the classical Runge-Kutta scheme of 4th order.

```cpp
  ponio::solver( lorenz_rhs, ponio::runge_kutta::rk_44(), u_ini, t_span, dt, obs);
```

## For more information

* [Documentation](https://ponio.readthedocs.io/en/latest/index.html)
* [Github repository](https://github.com/hpc-maths/ponio)
* [Demos](https://github.com/hpc-maths/ponio/tree/main/ponio/demos)
* [Notebooks examples](https://github.com/hpc-maths/ponio/tree/main/ponio/notebooks)
* [List of methods and their analysis](http://jmassot.perso.math.cnrs.fr/ponio/) (personal webpage of main developer)

## How to contribute

## Roadmap

## How to cite

## License

This project is licensed under the **BSD license**.

See [LICENSE](https://github.com/hpc-maths/ponio/blob/main/LICENSE) for more information.
