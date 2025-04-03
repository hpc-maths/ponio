# SciPy

> Fundamental algorithms for scientific computing in Python.

| Link           | [Web site](https://scipy.org/)                                                   |
|----------------|----------------------------------------------------------------------------------|
| Language       | Python                                                                           |
| License        | BSD-3-Clause license                                                             |
| Documentation  | [documentation](https://docs.scipy.org/doc/scipy/tutorial/index.html#user-guide) |
| Git repository | [@scipy](https://github.com/scipy/scipy)                                         |

## Lorenz equations

To solve Lorenz equations with a classical 4th order Runge-Kutta RK(4,4) we need to define it:

```{literalinclude} lorenz.py
  :lines: 6-17
```

and now call it with `scipy.integrate.solve_ivp` function:

```python
sol = solve_ivp(lorenz_system, t_span, y0, method=RK44, first_step=dt)
```

## Transport equation

SciPy provides mainly adaptive time step methods, so there is no classical 4th order Runge-Kutta RK(4,4) nor explicit Euler (or forward Euler or RK(1, 1)), we need to define it:

```{literalinclude} transport.py
  :lines: 6-13
```

and now call it with `scipy.integrate.solve_ivp` function:

```python
sol = solve_ivp(upwind, t_span, y0, method=Euler, first_step=dt)
```
