# ponio

> Presentation Overview of Numerical Integrator for ODE

La bibliothèque ponio est une collection d'intégrateurs en temps pour la résolution d'ODE et de PDE écrit en C++. ponio a vocation à être le plus simple d'utilisation sans fair de compromis sur les performances.

## Get started

### Lorenz equations

Intéressons nous à la résolution des équations de Lorenz :

$$
  \begin{cases}
    \dot{u}_1 &= \sigma(u_2 - u_1) \\
    \dot{u}_2 &= \rho u_1 - u_2 - u_1 u_3 \\
    \dot{u}_3 &= u_1 u_2 - \beta u_3
  \end{cases}
$$

ponio a été écrit pour résoudre des problèmes de la forme :

$$˜
  \dot{u} = f(t, u)
$$

The following steps describe how to solve this problem with samurai. It is important to note that these steps are generally the same whatever the equations we want to solve.

* Define the type of your state, generally named `state_t`. For the sake of simplicity, we will use a `std::valarray` container from the STL

```cpp
  using state_t = std::valarray<double>;
```

> **Becareful:** type `state_t` needs to support arithmetic operations. That why we don't use a `std:vector` in this example.

* Define a function, a lambda or a functor that represent the problem

```cpp
    const double sigma = 10.;
    const double rho   = 28.;
    const double beta  = 8./3.;

    auto lorenz = [=]( double t, state_t const& u ) -> state_t
    {
      return {
        sigma * ( u[1] - u[0] ),
        rho * u[0] - u[1] - u[0] * u[2],
        u[0] * u[1] - beta * u[2]
      };
    };
```

* Define $\Delta t$ and `t_span`

Quelque soit le problème, la première chose à faire sera d'écrire une fonction prenant le temps et un état (de type `state_t`), et renvoyant un état.

> **Attention :** Le type `state_t` doit pouvoir supporter des opérations arithmétiques.

For the sake of simplicity, nous prendrons ici le type `std::valarray<double>` pour le type de notre état. Pour les équations de Lorenz, la fonction $f:t,u \mapsto f(t,u)$ s'écrit :

```cpp
  std::valarray<double> lorenz( double t, std::valarray<double> const& u )
  {
    const double sigma = 10.;
    const double rho   = 28.;
    const double beta  = 8./3.;

    return {
        sigma * ( u[1] - u[0] ),
        rho * u[0] - u[1] - u[0] * u[2],
        u[0] * u[1] - beta * u[2]
    };
  }
```


You first need to define

```cpp
  #include <valarray>

  #include <solver/observer.hpp>
  #include <solver/runge_kutta.hpp>
  #include <solver/solver.hpp>
  #include <solver/time_span.hpp>

  int
  main()
  {
    using namespace observer;
    using state_t = std::valarray<double>;
    const double sigma = 10.;
    const double rho   = 28.;
    const double beta  = 8./3.;

    auto lorenz = [=]( double t, state_t const& u ) -> state_t
    {
      return {
        sigma * ( u[1] - u[0] ),
        rho * u[0] - u[1] - u[0] * u[2],
        u[0] * u[1] - beta * u[2]
      };
    };

    state_t u_ini = {1., 1., 1.};

    ponio::time_span<double> const t_span = {0., 20.};

    double const dt = 0.01;

    ponio::solve( lorenz, ponio::runge_kutta::rk_44(), u_ini, t_span, dt, "lorenz_sol.txt"_fobs );
  }
```

## Installation

### From conda


```
  conda install ponio
```

### From source

```
  git clone https://...
  cd ponio
  ...
```

## For more information

* [Documentation](#)
* [Github repository](#)
* [Demos](#)
* [Notebooks examples](#)
* [List of methods and their analysis](#)

## How to contribute

## Roadmap

## How to cite
