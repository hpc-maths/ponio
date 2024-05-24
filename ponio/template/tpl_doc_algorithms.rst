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

When the matrix :math:`A` is lower triangular with a diagonal, the Runge-Kutta method is called diagonal implicit (or DIRK). You have to provide a Jacobian function that returns the Jacobian matrix in point :math:`(t, u)` (see :cpp:class:`ponio::implicit_problem`). You can also provide an operator base definition (see :cpp:class:`ponio::implicit_operator_problem`).

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

.. note::

  Matrix :math:`A` and vector :math:`b` could contain some functions defined by

  .. math::

    \varphi_\ell(z) = \frac{e^z - \sum_{k=0}^{\ell-1} \frac{1}{k!}z^k }{z^\ell}


  and we use the notations :\math:`\varphi_\ell = \varphi_\ell(\Delta t L)` and :math:`\varphi_{\ell,j} = \varphi_\ell(c_j \Delta t L)`.


{% for rk in list_exprk %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endfor %}


Explicit stabilized Runge-Kutta methods
---------------------------------------

Some problems, like heat equation, require methods stabilized on the negative real axis. The ponio library provides a Runge-Kutta Chebyshev method of order 2, ROCK2 method (of order 2) and ROCK4 method (of order 4).

Runge-Kutta Chebyshev method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm of RKC2 is the following:

.. math::

   \begin{aligned}
      y_0 &= f(t^n, y^n) \\
      y_1 &= y_0 + \tilde{\mu}_1\Delta t f(t^n, y^n) \\
      y_j &= (1 - \mu_j - \nu_j)y^n + \mu_j y_{j-1} + \nu_j y_{j-2} \\
          &~~~~~~  + \tilde{\mu}_j \Delta t f(t^n + c_j\Delta t, y_{j-1}) + \tilde{\gamma}_j\Delta t f(t^n, y^n), \quad j=2,\dots, s
   \end{aligned}

with coefficients given by

.. math::

  \tilde{\mu}_1 = b_1 \omega_1

.. math::

  \mu_j = \frac{2b_j}{b_{j-1}}\omega_0, \quad \nu_j = -\frac{b_j}{b_{j-2}}, \quad \tilde{\mu}_j = \frac{2b_j}{b_{j-1}}\omega_1

.. math::

  \tilde{\gamma}_j = -(1 - b_{j-1}T_{j-1}(\omega_0))

where

.. math::

  b_0 = b_2, \quad b_1 = \frac{1}{\omega_0}, \quad
  b_j = \frac{T_j''(\omega_0)}{(T_j'(\omega_0))^2},\ j=2,\dots, s


.. math::

  c_0 = 0, \quad c_1 = c_2, \quad c_j = \frac{T_s'(\omega_0)}{T_s''(\omega_0)}\frac{T_j''(\omega_0)}{T_j'(\omega_0)}, \quad c_s = 1

and

.. math::

  \omega_0 = 1 + \frac{\epsilon}{s^2},\quad \omega_1 = \frac{T_s'(\omega_0)}{T_s''(\omega_0)},\quad \epsilon \approx \frac{2}{13}

and where :math:`T_j(x)` is the Chebyshev polynomial.


.. doxygenclass:: ponio::runge_kutta::chebyshev::explicit_rkc2
   :project: ponio
   :members:


ROCK2 method
~~~~~~~~~~~~

The algorithm of ROCK2 method is the following:

.. math::

  \begin{aligned}
    y_0 &= y^n \\
    y_1 &= y^n + \Delta t \mu_1 f(y^n) \\
    y_j &= \Delta t \mu_j f(y_{j-1}) - \nu_j y_{j-1} - \kappa_j y_{j-2}, \quad j=2,\dots, s-2\\
    y_{s-1} &= y_{s-2} + \Delta t \sigma f(u_{s-2}) \\
    y_{s}^\star &= y_{s-1} + \Delta t \sigma f(y_{s-1}) \\
    y^{n+1} &= y_s^\star - \Delta t \sigma (1 - \frac{\tau}{\sigma})( f(y_{s-1}) - f(y_{s-2}) )
  \end{aligned}

where :math:`\mu_j`, :math:`\nu_j` and :math:`\kappa_j` coefficients coming from a minimization problem.

.. doxygenclass:: ponio::runge_kutta::rock::rock2_impl
   :project: ponio
   :members:

.. doxygenfunction:: ponio::runge_kutta::rock::rock2(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::rock::rock2()
  :project: ponio


ROCK4 method
~~~~~~~~~~~~

The algorithm of ROCK4 method is the following:

.. math::

  \begin{aligned}
    y_0 &= y^n \\
    y_1 &= y^n + \Delta t \mu_1 f(y^n) \\
    y_j &= \Delta t \mu_j f(y_{j-1}) - \nu_j y_{j-1} - \kappa_j y_{j-2}, \quad j=2,\dots, s-4 \\
    y_{s-3} &= y_{s-4} + a_{21} \Delta t f(u_{s-4}) \\
    y_{s-2} &= y_{s-4} + a_{31} \Delta t f(u_{s-4}) + a_{32} \Delta t f(y_{s-3}) \\
    y_{s-1} &= y_{s-4} + a_{41} \Delta t f(u_{s-4}) + a_{42} \Delta t f(y_{s-3}) + a_{43} \Delta t f(y_{s-2}) \\
    y^{n+1} &= y_{s-4} + b_1 \Delta t f(y_{s-4}) + b_2 \Delta t f(y_{s-3}) + b_3 \Delta t f(y_{s-2}) + b_4 \Delta t f(y_{s-1})
  \end{aligned}

where :math:`\mu_j`, :math:`\nu_j` and :math:`\kappa_j` coefficients coming from a minimization problem, and :math:`a_{ij}`, :math:`b_i` coming from an order 4 method build with :math:`y_{s-4}` as initial condition.

.. doxygenclass:: ponio::runge_kutta::rock::rock4_impl
   :project: ponio
   :members:

.. doxygenfunction:: ponio::runge_kutta::rock::rock4(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::rock::rock4()
  :project: ponio
