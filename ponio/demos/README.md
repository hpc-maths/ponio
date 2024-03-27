# demos

The following table gives an overview over all examples.

| Section                                                                                                   | Brief Description                                                  | File                                                               |
|-----------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------|--------------------------------------------------------------------|
| [1. Arenstorf orbit](#1-arenstorf-orbit)                                                                  | This example shows how to use multiple adaptive time step methods  | [`arenstorf.cpp`](arenstorf.cpp)                                   |
| [2. Brownian movement](#2-brownian-movement)                                                              | This example shows how to use random in an ODE                     | [`brownian.cpp`](brownian.cpp)                                     |
| [3. Brusselator equations](#3-brusselator-equations)                                                      | The chemistry example of Brusselator (2 equations model)           | [`brusselator.cpp`](brusselator.cpp)                               |
| [4. Brusselator equations with DIRK method](#4-brusselator-equations-with-dirk-method)                    | This example shows hos to use DIRK methods                         | [`brusselator_dirk.cpp`](brusselator_dirk.cpp)                     |
| [5. Curtiss-Hirschfelder equation](#5-curtiss-hirschfelder-equation)                                      | This example shows how to use range and iterators on solution      | [`curtiss_hirschfelder.cpp`](curtiss_hirschfelder.cpp)             |
| [6. Curtiss-Hirschfelder equation with expRK method](#6-curtiss-hirschfelder-equation-with-exprk-method)  | This example shows how to use exponential Runge-Kutta methods      | [`curtiss_hirschfelder_exprk.cpp`](curtiss_hirschfelder_exprk.cpp) |
| [7. Exponential function](#7-exponential-function)                                                        | This example is the simplest example                               | [`exp.cpp`](exp.cpp)                                               |
| [8. Heat model](#8-heat-model)                                                                            | The classical heat equation solving with RKC2 method               | [`heat.cpp`](heat.cpp)                                             |
| [9. ROCK method](#9-rock-method)                                                                          | This example shows how to use ROCK2 and ROCK4 methods              | [`heat_rock.cpp`](heat_rock.cpp)                                   |
| [10. Samurai is hot](#10-samurai-is-hot)                                                                  | This example shows how to coupling ponio and samurai               | [`heat_samurai.cpp`](heat_samurai.cpp)                             |
| [11. Lorenz equations](#11-lorenz-equations)                                                              | The chaotic system example of Lorenz equations                     | [`lorenz.cpp`](lorenz.cpp)                                         |
| [12. Lorenz equations with multiple methods](#12-lorenz-equations-with-multiple-methods)                  | This example shows how to use splitting methods and Lawson methods | [`lorenz_tuto.cpp`](lorenz_tuto.cpp)                               |
| [13. Lotka-Volterra model](#13-lotka-volterra-model)                                                      | The classical predator–prey model of Lotka-Volterra                | [`lotka_volterra.cpp`](lotka_volterra.cpp)                         |
| [14. Nagumo equation](#14-nagumo-equation)                                                                | This example shows how to use ponio to mesure order of a method    | [`nagumo.cpp`](nagumo.cpp)                                         |
| [15. Pendulum equation](#15-pendulum-equation)                                                            | The classical pendulum equation                                    | [`pendulum.cpp`](pendulum.cpp)                                     |

## 1. Arenstorf orbit

The system of differential equations for the Arenstorf orbit are:

$$
  \begin{cases}
    \ddot{x} &= x + 2\dot{y} - \frac{1-\mu}{r_1^3}(x+\mu) - \frac{\mu}{r_2^3}(x-1+\mu) \\
    \ddot{y} &= y - 2\dot{x} - \frac{1-\mu}{r_1^3}y - \frac{\mu}{r_2^3}y
  \end{cases}
$$

where

$$
  r_1 = \sqrt{(x+\mu)^2 + y^2},\quad r_2 = \sqrt{(x-1+\mu)^2 + y^2}
$$

parameter $\mu=0.012277471$ and the initial condition gives by:

$$
  x(0) = 0.994,\quad \dot{x}(0) = 0,\quad y(0) = 0,\quad \dot{y}(0) = -2.00158510637908252240537862224
$$

To solve this kind of problem with ponio, first of all you should rewrite it as the form: $\dot{u} = f(t, u)$, here we classically take

$$
  u = \begin{pmatrix}
    x \\
    y \\
    \dot{x} \\
    \dot{y}
  \end{pmatrix}
$$

So we have:

$$
  \dot{u} = \begin{pmatrix}
    \dot{x} \\
    \dot{y} \\
    \ddot{x} \\
    \ddot{y}
  \end{pmatrix} = \begin{pmatrix}
    \dot{x} \\
    \dot{y} \\
    x + 2\dot{y} - \frac{1-\mu}{r_1^3}(x+\mu) - \frac{\mu}{r_2^3}(x-1+\mu) \\
    y - 2\dot{x} - \frac{1-\mu}{r_1^3}y - \frac{\mu}{r_2^3}y
  \end{pmatrix} = f(t, u)
$$

In this example we solve this system with some explicit adaptive time step methods from [Dormand, J.R., Prince, P.J., A family of embedded Runge-Kutta formulae (1980) *Journal of Computational and Applied Mathematics*](http://dx.doi.org/10.1016/0771-050x(80)90013-3)

| Arenstorf orbit                                  | Arenstorf velocity                                  |
|--------------------------------------------------|-----------------------------------------------------|
| ![Arenstorf orbit](img/1-arenstorf-orbit_01.png) | ![Arenstorf velocity](img/1-arenstorf-orbit_02.png) |

## 2. Brownian movement

$$
  \begin{cases}
    \dot{x} = X(t) \\
    \dot{y} = Y(t)
  \end{cases}
$$

## 3. Brusselator equations

$$
  \begin{cases}
    \dot{x} &= m_a - (m_b + 1)x + x^2y \\
    \dot{y} &= m_b x - x^2y
  \end{cases}
$$

## 4. Brusselator equations with DIRK method

$$
  \begin{cases}
    \dot{x} &= m_a - (m_b + 1)x + x^2y \\
    \dot{y} &= m_b x - x^2y
  \end{cases}
$$

## 5. Curtiss-Hirschfelder equation

$$
  \dot{y} = k(\cos(t) - y)
$$

## 6. Curtiss-Hirschfelder equation with expRK method

$$
  \dot{y} = k(\cos(t) - y)
$$

## 7. Exponential function

$$
  \dot{y} = y
$$

## 8. Heat model

$$
  \dot{u} = -\partial_{xx} u
$$

## 9. ROCK method

$$
  \dot{u} = -\partial_{xx} u
$$

## 10. Samurai is hot

$$
  \dot{u} = -\partial_{xx} u
$$

## 11. Lorenz equations

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - xz \\
    \dot{z} &= xy - \beta z
  \end{cases}
$$

## 12. Lorenz equations with multiple methods

$$
  \begin{cases}
    \dot{x} &= \sigma (y - x) \\
    \dot{y} &= \rho x - y - xz \\
    \dot{z} &= xy - \beta z
  \end{cases}
$$

## 13. Lotka-Volterra model

$$
  \begin{cases}
    \dot{x} = \alpha x - \beta xy \\
    \dot{y} = \delta xy - \gamma y
  \end{cases}
$$

## 14. Nagumo equation

$$
  \partial_t u = d \partial_{xx}u + ku^2(1-u)
$$

## 15. Pendulum equation

$$
  \ddot{\theta} + b\dot{\theta} + c\sin(\theta) = 0
$$
