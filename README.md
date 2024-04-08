# ponio

> Presentation Overview of Numerical Integrator for ODE

The library ponio is a collection of time integrators for solving differential equations written in C++. The purpose of this library is to supply efficient and flexible C++ implementation of solvers for various differential equations. Main method classes are :

* explicit Runge-Kutta methods (eRK)
* diagonal implicit Runge-Kutta methods (DIRK)
* Lawson methods (based with an underlying Runge-Kutta method) (LRK)
* exponential Runge-Kutta methods (expRK)
* Runge-Kutta Chebyshev (RKC)
* splitting method (Lie or Strang)

This library aims to be the easiest to use without compromising on performance.

<details>
<summary>Table of Contents</summary>

- [Get started](#get-started)
- [Features](#features)
- [Installation](#installation)
  - [From conda](#from-conda)
  - [From source](#from-source)
- [For more information](#for-more-information)
- [How to contribute](#how-to-contribute)
- [License](#license)

</details>

## Get started

### Lorenz equations

In this section, we will present how to solve the Lorenz equations, defined:

$$
  \begin{cases}
    \dot{y}_1 &= \sigma(y_2 - y_1) \\
    \dot{y}_2 &= \rho y_1 - y_2 - y_1 y_3 \\
    \dot{y}_3 &= y_1 y_2 - \beta y_3
  \end{cases}
$$

This model is weell know to be a chaotic system, here the solution plots form every 92 methods provide by ponio.

| Lorenz attractor |
|------------------|
| ![Lorenz attractor](https://private-user-images.githubusercontent.com/7198360/318870562-f263304a-6d1e-4698-9be7-963bad655cb8.gif?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTI1NTIxMTQsIm5iZiI6MTcxMjU1MTgxNCwicGF0aCI6Ii83MTk4MzYwLzMxODg3MDU2Mi1mMjYzMzA0YS02ZDFlLTQ2OTgtOWJlNy05NjNiYWQ2NTVjYjguZ2lmP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQVZDT0RZTFNBNTNQUUs0WkElMkYyMDI0MDQwOCUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyNDA0MDhUMDQ1MDE0WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9MjcwZWMwNTIwNTY0ZmVmZWMzMzA4ZWZmNmY2MTE3ZGUwNDIzNjJkN2M5ZTY5ODQ3M2EyNDE4NGY3OWIwOTUxYiZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QmYWN0b3JfaWQ9MCZrZXlfaWQ9MCZyZXBvX2lkPTAifQ.anvhY6CgxFFdmoScWrHBxQeYHnecVbWaTc_P5HSRPnU) |


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

  auto obs = "sol.txt"_fobs;
```

This line will create a file `sol.txt` and push all output data in it. ponio observes write data as `csv` format : `current_time, current_state, current_time_step`, with eventually multiple columns for `current_state`.

Now you are ready to solve your problem with an explicit Runge-Kutta method, see [algorithm overview](https://ponio.readthedocs.io/en/latest/api/algorithm.html). In this example we will use the classical Runge-Kutta scheme of 4th order.

```cpp
  ponio::solver( lorenz_rhs, ponio::runge_kutta::rk_44(), u_ini, t_span, dt, obs);
```

The whole example can be found in [here](https://github.com/hpc-maths/ponio/blob/main/ponio/examples/lorenz.cpp).


### More examples

More examples can be found in [notebooks](https://github.com/hpc-maths/ponio/tree/main/ponio/notebooks) or [examples](https://github.com/hpc-maths/ponio/tree/main/ponio/examples) directories.

## Features

* [x] Explicit Runge-Kutta methods from their Butcher tableau
* [x] Diagonal-implicit Runge-Kutta methods from their Butcher tableau
* [x] Lawson methods from all explicit Runge-Kutta methods
* [x] exponential Runge-Kutta methods from their Butcher tableau
* [x] Runge-Kutta Chebyshev method of order 2
* [x] ROCK 2 and ROCK 4 methods
* [x] Splitting methods : Lie splitting method and Strang splitting method
* [x] PIROCK method
* [ ] Additive Runge-Kutta methods (IMEX) from their Butcher tableau
* [x] Coupling ponio and adaptive mesh library [samurai](https://github.com/hpc-maths/samurai)
* [x] Coupling ponio and linear algebra library [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* [ ] Parareal method
* [ ] Simplify multi-equations problem

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
```

* With [pixi](https://pixi.sh/latest/)

This method will install all dependencies for all examples.

```
  pixi install
  pixi build
```

* With only cmake

Get only sources to run a project

```
  cmake . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build --target install
```

## For more information

* [Documentation](https://ponio.readthedocs.io/en/latest/index.html)
* [Github repository](https://github.com/hpc-maths/ponio)
* [Examples](https://github.com/hpc-maths/ponio/tree/main/ponio/examples)
* [Notebooks examples](https://github.com/hpc-maths/ponio/tree/main/ponio/notebooks)
* [List of methods and their analysis](http://jmassot.perso.math.cnrs.fr/ponio/) (personal webpage of main developer)

## How to contribute

First off, thanks for taking the time to contribute! Contributions are what make the open-source community such an amazing place to learn, inspire, and create. Any contributions you make will benefit everybody else and are greatly appreciated.

Please read [our contribution guidelines](https://github.com/hpc-maths/ponio/blob/main/ponio/doc/CONTRIBUTING.md), and thank you for being involved!

## License

This project is licensed under the **BSD license**.

See [LICENSE](https://github.com/hpc-maths/ponio/blob/main/LICENSE) for more information.
