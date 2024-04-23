Algorithms
==========

List of all algorithms (skeleton of each numerical method to solve)

Runge-Kutta methods
-------------------

A Runge-Kutta method is defined in ponio by its Butcher tableau

.. math::

   \begin{array}{c|c}
      c & A \\
      \hline
        & b^\top
   \end{array}

which is stored as a JSON file in ``database`` folder.

Runge-Kutta methods are useful to solve a problem like

.. math::

   \dot{u} = f(t, u)

The method, with the previous Butcher tableau, reads as

.. math::

   \begin{aligned}
      u^{(i)} &= u^n + \Delta t \sum_j a_{ij} k_j, \quad i = 1, \dots, s \\
      k_i &= f(t^n + c_i\Delta t, u^{(i)}) \\
      u^{n+1} &= u^n + \Delta t \sum_i b_i k_i
   \end{aligned}


Explicit methods
~~~~~~~~~~~~~~~~

When the matrix :math:`A` is strictly lower triangular, the Runge-Kutta method is called explicit. The only think you need to provide to solve a problem with this kind of method is the function :math:`f`. See below for the list of explicit Runge-Kutta methods in ponio.

{% for rk in list_erk %}{% if rk.b2 is undefined %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

The ponio library provides also adaptive time step methods.

{% for rk in list_erk %}{% if rk.b2 is defined %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Diagonal implicit methods
~~~~~~~~~~~~~~~~~~~~~~~~~

When the matrix :math:`A` is lower triangular with a diagonal, the Runge-Kutta method is called diagonal implicit (or DIRK). You have to provide a Jacobian function that returns the Jacobian matrix in point :math:`(t, u)` (see :ref:`implcit_problem`). You can also provide an operator base definition (see :ref:`implicit_operator_problem`).

{% for rk in list_dirk %}
.. doxygenfunction:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endfor %}


Lawson methods
--------------

Lawson methods are build to solve a problem with a linear and nonlinear part, to solve exactly the problem when the nonlinear part goes to zero. This class of problem can be write as

.. math::

   \dot{u} = L u + N(t, u)

First, we introduce the change of variable

.. math::

   v(t) = e^{-Lt}u(t)

which yields the equation

.. math::

   \dot{v}(t) = -Le^{-Lt}u(t) + e^{-Lt}\dot{u}(t)

which can be rewrite in therm of :math:`v` as

.. math::

   \dot{v}(t) = e^{-Lt}N(t, e^{-Lt}v) = \tilde{N}(t, v)

We solve this equation with a classical Runge-Kutta method RK(:math:`s`, :math:`n`) with :math:`s` stages and of order :math:`n`. We rewrite the scheme in terms of the variable :math:`u`, which yields the following scheme

.. math::

   \begin{aligned}
      u^{(i)} &= u^n + \Delta t \sum_j a_{ij} k_j, \quad i = 1, \dots, s \\
      k_i &= e^{-c_i\Delta t L} N(t^n + c_i\Delta t, e^{c_i\Delta t L} u^{(i)}) \\
      u^{n+1} &= e^{\Delta t L}\left( u^n + \Delta t \sum_i b_i k_i \right)
   \end{aligned}


In ponio, Lawson methods have the same name of the underlying Runge-Kutta method prefixed by ``l``.


Explicit methods
~~~~~~~~~~~~~~~~

{% for rk in list_erk %}{% if rk.b2 is undefined %}
.. doxygenvariable:: ponio::runge_kutta::l{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

{% for rk in list_erk %}{% if rk.b2 is defined %}
.. doxygenvariable:: ponio::runge_kutta::l{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Exponential Runge-Kutta methods
-------------------------------

Like Lawson methods, exponential Runge-Kutta methods are build to solve a problem with a linear and non-linear part, to solve exactly the problem when the nonlinear part goes to zero. This class of problem can be write as

.. math::

   \dot{u} = L u + N(t, u)

We solve this on a time step between :math:`0` and :math:`\Delta t`

.. math::

  u(t^n + \Delta t) = e^{\Delta t L} U + \int_0^{\Delta t} e^{(\Delta t - s)L}N(t^n + s, u(t^n+s)) \,\mathrm{d}s

Interpolation of the integral yields to build a custom Runge-Kutta method which reads as

.. math::

    \begin{aligned}
        u^{(i)} &= u^n + \Delta t\sum_j a_{ij}(\Delta t L)\cdot ( k_j + Lu^n ), \quad i = 1, \dots, s \\
        k_i     &= N(t^n + c_i\Delta t, u^{(i)}) \\
        u^{n+1} &= u^n + \Delta t \sum_i b_i(\Delta t L)\cdot( k_i + Lu^n)
    \end{aligned}


{% for rk in list_exprk %}
.. doxygenvariable:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endfor %}


Runge-Kutta Chebyshev methods
-----------------------------

.. doxygenclass:: ponio::runge_kutta::chebyshev::explicit_rkc2
   :project: ponio
   :members:
