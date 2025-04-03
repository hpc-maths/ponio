# GSL - GNU Scientific Library

> The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free software under the GNU General Public License.

| Link           | [Web site](https://www.gnu.org/software/gsl/)              |
|----------------|------------------------------------------------------------|
| Language       | C                                                          |
| License        | GPL                                                        |
| Documentation  | [documentation](https://www.gnu.org/software/gsl/doc/html) |

## Transport equation

To define the explicit Euler method, we define a ``euler_state_t`` data structure to store current state and intermediate step:

```{literalinclude} transport.c
  :lines: 14-19
```

and add a new ``gsl_odeiv2_step_type`` with following lines:

```{literalinclude} transport.c
  :lines: 157-167
```

where we should define how to allocate a ``euler_state_t`` (``euler_alloc``), increment it with ``euler_apply``, reset data with ``euler_reset``, a function to return the order ``euler_order`` and a deallocator ``euler_free``.
