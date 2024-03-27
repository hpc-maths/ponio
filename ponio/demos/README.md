# demos

The following table gives an overview over all examples.

| File                                                               | Brief Description                                                  | Section                                                                                                   |
|--------------------------------------------------------------------|--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| [`arenstorf.cpp`](arenstorf.cpp)                                   | This example shows how to use multiple adaptive time step methods  | [1. Arenstorf orbit](#1-arenstorf-orbit)                                                                  |
| [`brownian.cpp`](brownian.cpp)                                     | This example shows how to use random in an ODE                     | [2. Brownian movement](#2-brownian-movement)                                                              |
| [`brusselator.cpp`](brusselator.cpp)                               | The chemistry example of Brusselator (2 equations model)           | [3. Brusselator equations](#3-brusselator-equations)                                                      |
| [`brusselator_dirk.cpp`](brusselator_dirk.cpp)                     | This example shows hos to use DIRK methods                         | [4. Brusselator equations with DIRK method](#4-brusselator-equations-with-dirk-method)                    |
| [`curtiss_hirschfelder.cpp`](curtiss_hirschfelder.cpp)             | This example shows how to use range and iterators on solution      | [5. Curtiss-Hirschfelder equation](#5-curtiss-hirschfelder-equation)                                      |
| [`curtiss_hirschfelder_exprk.cpp`](curtiss_hirschfelder_exprk.cpp) | This example shows how to use exponential Runge-Kutta methods      | [6. Curtiss-Hirschfelder equation with expRK method](#6-curtiss-hirschfelder-equation-with-exprk-method)  |
| [`exp.cpp`](exp.cpp)                                               | This example is the simplest example                               | [7. Exponential function](#7-exponential-function)                                                        |
| [`heat.cpp`](heat.cpp)                                             | The classical heat equation solving with RKC2 method               | [8. Heat model](#8-heat-model)                                                                            |
| [`heat_rock.cpp`](heat_rock.cpp)                                   | This example shows how to use ROCK2 and ROCK4 methods              | [9. ROCK method](#9-rock-method)                                                                          |
| [`heat_samurai.cpp`](heat_samurai.cpp)                             | This example shows how to coupling ponio and samurai               | [10. Samurai is hot](#10-samurai-is-hot)                                                                  |
| [`lorenz.cpp`](lorenz.cpp)                                         | The chaotic system example of Lorenz equations                     | [11. Lorenz equations](#11-lorenz-equations)                                                              |
| [`lorenz_tuto.cpp`](lorenz_tuto.cpp)                               | This example shows how to use splitting methods and Lawson methods | [12. Lorenz equations with multiple methods](#12-lorenz-equations-whith-multiple-methods)                 |
| [`lotka_volterra.cpp`](lotka_volterra.cpp)                         | The classical predator–prey model of Lotka-Volterra                | [13. Lotka-Volterra model](#13-lotka-voleterra-model)                                                     |
| [`nagumo.cpp`](nagumo.cpp)                                         | This example shows how to use ponio to mesure order of a method    | [14. Nagumo equation](#14-nagumo-model)                                                                   |
| [`pendulum.cpp`](pendulum.cpp)                                     | The classical pendulum equation                                    | [15. Pendulum equation](#15-pendulum-equation)                                                            |

## 1. Arenstorf orbit

$$
  \begin{cases}
    \ddot{x} &= x + 2\dot{y} - \frac{1-\mu}{r_1^3}(x+\mu) - \frac{\mu}{r_2^3}(x-1+\mu) \\
    \ddot{y} &= y - 2\dot{x} - \frac{1-\mu}{r_1^3}y - \frac{\mu}{r_2^3}y
  \end{cases}
$$

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
  \dot{u} = \partial_{xx} u
$$

## 9. ROCK method

$$
  \dot{u} = \partial_{xx} u
$$

## 10. Samurai is hot

$$
  \dot{u} = \partial_{xx} u
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
  \ddot{\theta} + b\dot{theta} + c\sin(\theta) = 0
$$
