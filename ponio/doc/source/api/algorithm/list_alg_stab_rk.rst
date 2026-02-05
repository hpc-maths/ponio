List of extended stabilized Runge-Kutta methods
===============================================

Some problems, like heat equation, require methods stabilized on the negative real axis. The ponio library provides a Runge-Kutta Chebyshev method of order 2, ROCK2 method (of order 2), ROCK4 method (of order 4) and a Runge-Kutta Legendre method of order 1 or 2.


Runge-Kutta Chebyshev method
----------------------------

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

.. warning::

  In ponio library the number of stages of RKC2 is static, you can't get dynamic number of stages for this method, look at :cpp:func:`ponio::runge_kutta::rock::rock2` or :cpp:func:`ponio::runge_kutta::rock::rock4` methods for dynamic number of stages (and optional adaptive time step method).


ROCK2 method
------------

We write the method ROCK2 presented in :cite:`abdulle:2001`. The algorithm of ROCK2 method is the following:

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

Helper functions
~~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::runge_kutta::rock::rock2(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::rock::rock2()
  :project: ponio


ROCK4 method
------------

We write the method ROCK2 presented in :cite:`abdulle:2002`. The algorithm of ROCK4 method is the following:

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

Helper functions
~~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::runge_kutta::rock::rock4(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::rock::rock4()
  :project: ponio


Runge-Kutta Legendre method
---------------------------

An other way to get a stabilized Runge-Kutta method is to use Legendre polynomials, we follow presentation in :cite:`meyer:2014`.

Runge-Kutta Legend first order method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The algorithm of RKL1 is the following:

.. math::

   \begin{aligned}
      y^{(0)} &= y^n \\
      y^{(1)} &= y^n + \tilde{\mu}_1\Delta t f(t^n, y^n) \\
      y^{(j)} &= \mu_jy^{(j-1)} + \nu_j y^{(j-2)} + \tilde{\mu}_j\Delta t f(t^n, y^{(j-1)}), \quad j=2,\dots, s \\
      y^{(n+1)} &= y^{(s)}
   \end{aligned}

with coefficients given by

.. math::

  \mu_j = \frac{2j-1}{j}, \qquad \nu_j = \frac{1-j}{j}

.. math::

  \tilde{\mu}_j = \frac{2j-1}{j}\frac{2}{s^2 + s}

where :math:`s` is the number of stages of the method.

.. doxygenclass:: ponio::runge_kutta::legendre::explicit_rkl1
   :project: ponio
   :members:

.. warning::

  In ponio library the number of stages of RKL1 is static, you can't get dynamic number of stages for this method, look at :cpp:func:`ponio::runge_kutta::rock::rock2` or :cpp:func:`ponio::runge_kutta::rock::rock4` methods for dynamic number of stages (and optional adaptive time step method).


Runge-Kutta Legend second order method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm of RKL2 is the following:

.. math::

   \begin{aligned}
      y^{(0)} &= y^n \\
      y^{(1)} &= y^{(0)} + \tilde{\mu}_1\Delta t f(t^n, y^{(0)}) \\
      y^{(j)} &= \mu_jy^{(j-1)} + \nu_j y^{(j-2)} + (1-\mu_j-\nu_j)y^{(0)} + \tilde{\mu}_j\Delta t f(t^n, y^{(j-1)}) + \tilde{\gamma}_j\Delta tf(t^n, y^{(0)}), \quad j=2,\dots, s \\
      y^{(n+1)} &= y^{(s)}
   \end{aligned}

with coefficients given by

.. math::

  \mu_j = \frac{2j-1}{j}\frac{b_j}{b_{j-1}}, \qquad \nu_j = -\frac{j-1}{j}\frac{b_j}{b_{j-2}}

.. math::

  \tilde{\mu}_j = \mu_j w_1,\ 1<j \qquad  \tilde{\mu}_1 = b_1 w_1

.. math::

  \tilde{\gamma}_j = -a_{j-1}\tilde{\mu}_j

where

.. math::

  b_0 = b_1 = b_2 = \frac{1}{3},\qquad b_j = \frac{j^2 + j - 2}{2j(j+1)},\ 2\leq j

and

.. math::

  a_j = 1 - b_j, \qquad w_1 = \frac{4}{s^2 + s - 2}


where :math:`s` is the number of stages of the method

.. doxygenclass:: ponio::runge_kutta::legendre::explicit_rkl2
   :project: ponio
   :members:

.. warning::

  In ponio library the number of stages of RKL2 is static, you can't get dynamic number of stages for this method, look at :cpp:func:`ponio::runge_kutta::rock::rock2` or :cpp:func:`ponio::runge_kutta::rock::rock4` methods for dynamic number of stages (and optional adaptive time step method).
