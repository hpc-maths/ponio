Time integrators
================

The ponio library provides different category of numerical integrators:

* Runge-Kutta methods based on Butcher tableaus,
* Extended stability methods,
* Splitting methods,
* IMEX methods with extended stability method.

For all examples in this section we will solve the Curtiss-Hirschfelder equation

.. math::

  \dot{y} = k(\cos(t) - y)

with stiff parameter :math:`k=50`, initial condition :math:`y(0) = y_0 = 2` and time :math:`t\in[0, 4]` with an initial time step :math:`\Delta t=0.05`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 30-34
  :lineno-start: 30
  :linenos:

Runge-Kutta methods from Butcher tableau
----------------------------------------

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

When the matrix :math:`A` is strictly lower triangular, the Runge-Kutta method is called explicit. The only think you need to provide to solve a problem with this kind of method is the function :math:`f`.

.. seealso::

   See the :doc:`list of explicit Runge-Kutta methods <../api/algorithm/list_alg_erk>` in ponio.

For an explicit method, we can only define a the function :math:`f` as a lambda as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 37-40
  :lineno-start: 37
  :linenos:

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 42
  :lineno-start: 42
  :linenos:

Embedded methods
~~~~~~~~~~~~~~~~

Some Butcher tableaus provide a second :math:`b` line:

.. math::

   \begin{array}{c|c}
      c & A \\
      \hline
        & b^\top
        & \hat{b}^\top
   \end{array}

The Butcher tableau reads as

.. math::

   \begin{aligned}
      u^{(i)} &= u^n + \Delta t \sum_j a_{ij} k_j, \quad i = 1, \dots, s \\
      k_i &= f(t^n + c_i\Delta t, u^{(i)}) \\
      u^{n+1} &= u^n + \Delta t \sum_i b_i k_i \\
      \tilde{u}^{n+1} &= u^n + \Delta t \sum_i \hat{b}_i k_i
   \end{aligned}

with :math:`u^{n+1}` and :math:`\tilde{u}^{n+1}` two estimate of solution at time :math:`t^{n+1}` of different order, we can compute a normalized error estimate with the following formula, with a given absolute tolerance :math:`a_{tol}` respectively relative tolerance :math:`r_{tol}`

.. math::

  e = \sqrt{ \frac{1}{N}\sum_i \left( \frac{|u^{n+1}_i - \tilde{u}^{n+1}_i|}{a_{tol} + r_{tol} \max( |u^n_i|, |u^{n+1}_i| )} \right)^2 }

This function is defined in :cpp:func:`ponio::detail::error_estimate`. Now compare the normalized error to :math:`1` to accept or not the step, and compute a new time step

.. math::

  \begin{aligned}
    \Delta t_\text{new} &= 0.9 \sqrt[n]{ \frac{tol}{e} } \Delta t \\
    \Delta t            &= \min\left( \max\left( \Delta t_\text{new}, 0.2\Delta t \right) , 5\Delta t \right)
  \end{aligned}

Most common embedded Runge-Kutta methods come from :cite:`dormand:1980` and :cite:`prince:1981`, that why tey are sometime call Dormand-Prince methods.

.. seealso::

   See the :doc:`list of embedded Runge-Kutta methods <../api/algorithm/list_alg_dp>` in ponio.

The most common embedded methods are explicit methods, so we can only define a the function :math:`f` as a lambda as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 46-49
  :lineno-start: 46
  :linenos:

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 51
  :lineno-start: 51
  :linenos:

Diagonal implicit methods
~~~~~~~~~~~~~~~~~~~~~~~~~

When the matrix :math:`A` is lower triangular with a diagonal, the Runge-Kutta method is called diagonal implicit (or DIRK). You have to provide a Jacobian function that returns the Jacobian matrix in point :math:`(t, u)` (see :cpp:class:`ponio::implicit_problem`). You can also provide an operator base definition (see :cpp:class:`ponio::implicit_operator_problem`).

.. seealso::

   See the :doc:`list of diagonal implicit Runge-Kutta methods <../api/algorithm/list_alg_dirk>` in ponio.

For an diagonal implicit Runge-Kutta method, we need the function :math:`f` and also its Jacobian function :math:`\mathrm{d}f`, and store them into a :cpp:class:`ponio::implicit_problem`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 55-64
  :lineno-start: 55
  :linenos:

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 66
  :lineno-start: 66
  :linenos:


Lawson methods
~~~~~~~~~~~~~~

Lawson methods was initially propose into :cite:`lawson:1967`. Lawson methods are build to solve a problem with a linear and nonlinear part, to solve exactly the problem when the nonlinear part goes to zero. This class of problem can be write as

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


In ponio, Lawson methods have the same name of the underlying explicit Runge-Kutta method prefixed by ``l``.

.. seealso::

   See the :doc:`list of Lawson Runge-Kutta methods <../api/algorithm/list_alg_lrk>` in ponio.

For a problem split into a linear and nonlinear part, we need to define linear part as a scalar or a matrix and the nonlinear part as a function, and store them into a :cpp:class:`ponio::laxson_problem`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 70-76
  :lineno-start: 70
  :linenos:

Because sometime we want to define exponential function with a Pade approximant, we need to specify an exponential function. In this case we use :code:`std::exp<double>` standard function. And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 78
  :lineno-start: 78
  :linenos:

Exponential Runge-Kutta methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


  and we use the notations :math:`\varphi_\ell = \varphi_\ell(\Delta t L)` and :math:`\varphi_{\ell,j} = \varphi_\ell(c_j \Delta t L)`.

.. seealso::

   See the :doc:`list of exponential Runge-Kutta methods <../api/algorithm/list_alg_exprk>` in ponio.

For a problem split into a linear and nonlinear part, we need to define linear part as a scalar or a matrix and the nonlinear part as a function, and store them into a :cpp:class:`ponio::laxson_problem`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 82-88
  :lineno-start: 82
  :linenos:

 And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 90
  :lineno-start: 90
  :linenos:

----


Extended stability methods
--------------------------

Some problems, like heat equation, require methods stabilized on the negative real axis. The ponio library provides a Runge-Kutta Chebyshev method of order 2, ROCK2 method (of order 2), ROCK4 method (of order 4) and a Runge-Kutta Legendre method of order 1 or 2.

.. seealso::

   See the :doc:`list of extended stability Runge-Kutta methods <../api/algorithm/list_alg_stab_rk>` in ponio.


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

Runge-Kutta Chebyshev method is an explicit method, so we only need to define :math:`f` function as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 94-97
  :lineno-start: 94
  :linenos:

We need to specify how many stages we want (we choose 5), and next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 99
  :lineno-start: 99
  :linenos:


ROCK2 method
~~~~~~~~~~~~

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

ROCK2 method is an explicit method, so we only need to define :math:`f` function as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 103-106
  :lineno-start: 103
  :linenos:

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 108
  :lineno-start: 108
  :linenos:

.. hint::

   The ROCK2 method computes dynamically the number of stages with power method to estimate the spectral radius of your operator, but you can provide an other estimator of spectral radius if needed.

ROCK4 method
~~~~~~~~~~~~

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

ROCK4 method is an explicit method, so we only need to define :math:`f` function as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 112-115
  :lineno-start: 112
  :linenos:

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 117
  :lineno-start: 117
  :linenos:

.. hint::

   The ROCK4 method computes dynamically the number of stages with power method to estimate the spectral radius of your operator, but you can provide an other estimator of spectral radius if needed.


Runge-Kutta Legendre method
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An other way to get a stabilized Runge-Kutta method is to use Legendre polynomials, we follow presentation in :cite:`meyer:2014`.

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

Runge-Kutta Legendre first-order method is an explicit method, so we only need to define :math:`f` function as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 121-124
  :lineno-start: 121
  :linenos:

We need to specify how many stages we want (we choose 5), and next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 126
  :lineno-start: 126
  :linenos:

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

Runge-Kutta Legendre second-order method is an explicit method, so we only need to define :math:`f` function as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 130-133
  :lineno-start: 130
  :linenos:

We need to specify how many stages we want (we choose 5), and next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 135
  :lineno-start: 135
  :linenos:

----


Splitting methods
-----------------

Some methods to solve a problem of the form:

.. math::

  \dot{u} = \sum_i f_i(t,u)

with the initial condition :math:`u(t=0)=u_0`. We note with :math:`\phi_{\tau}^{[f_i]}(t^n,\tilde{u}^n)` the solution at time :math:`t^n+\tau` of the subproblem :math:`\dot{u}=f_i(t,u)` with the initial condition :math:`u(t^n)=\tilde{u}^n`.

.. seealso::

   See the :doc:`list of splitting methods <../api/algorithm/splitting>` in ponio.


Lie Splitting method
~~~~~~~~~~~~~~~~~~~~

In Lie splitting method, the solution is computed as:

.. math::

   u^{n+1} = \mathcal{L}^{\Delta t}(t^n, u^n) = \phi_{\Delta t}^{[f_1]}\circ \cdots \circ \phi_{\Delta t}^{[f_k]} (t^n,u^n)

For splitting methods in ponio, we first need to define all sub-problem (here only :math:`f_1` and :math:`f_2`) and store them into a :cpp:class:`ponio::problem`

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 139-148
  :lineno-start: 139
  :linenos:

We also define the Lie splitting method and the method to solve each substep. In this example we use two RK(3,3) methods, and we specify the time step for subcycle

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 150-151
  :lineno-start: 150
  :linenos:

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 153
  :lineno-start: 153
  :linenos:


Strang Splitting method
~~~~~~~~~~~~~~~~~~~~~~~

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \mathcal{S}^{\Delta t}(t^n, u^n) = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}
              \circ \phi_{\Delta t}^{[f_k]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

For splitting methods in ponio, we first need to define all sub-problem (here only :math:`f_1` and :math:`f_2`) and store them into a :cpp:class:`ponio::problem`

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 157-166
  :lineno-start: 157
  :linenos:

We also define the Strang splitting method and the method to solve each substep. In this example we use two RK(3,3) methods, and we specify the time step for subcycle

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 168-169
  :lineno-start: 168
  :linenos:

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder.cpp
  :language: cpp
  :lines: 171
  :lineno-start: 171
  :linenos:

Adaptive time step Strang splitting method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \mathcal{S}^{\Delta t}(t^n, u^n) = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}
              \circ \phi_{\Delta t}^{[f_k]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

The approximation of order 1 is computed with a shifted Strang splitting method:

.. math::

   \tilde{u}^{n+1} = \mathcal{S}_{\delta}^{\Delta t}(t^n, u^n) = \phi_{(\frac{1}{2}-\delta)\Delta t}^{[f_1]}\circ\phi_{\frac{\Delta t}{2}}^{[f_2]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}
              \circ \phi_{\Delta t}^{[f_k]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}\circ\cdots\circ\phi_{\frac{\Delta t}{2}}^{[f_2]}\circ \phi_{(\frac{1}{2}+\delta)\Delta t}^{[f_1]}
              (t^n,u^n)

with :math:`\delta\in[-1/2, 0)\cup(0,1/2]`

.. note::

   Only the first (and last) step is shifted by :math:`\delta\Delta t`.

The difficulty in adaptive time step Strang splitting method is to find a good shift coefficient :math:`\delta` for a given tolerance :math:`\eta`. The complete strategy of an iteration provides in :cite:`descombes:2015` is

1. With a given tolerance :math:`\eta`, time step :math:`\Delta t` and shift :math:`\delta` compute estimates

   .. math::

      \begin{aligned}
         u^{n+1} &= \mathcal{S}^{\Delta t}(t^n, u^n) \\
         \tilde{u}^{n+1} &= \mathcal{S}_{\delta}^{\Delta t}(t^n, u^n)
      \end{aligned}
2. Compute new time step from local error

   .. math::

      \Delta t^{new} = s \Delta t \sqrt{\frac{\eta}{\| u^{n+1} - \tilde{u}^{n+1} \|}}
3. If local error is lower than tolerance: :math:`\| u^{n+1} - \tilde{u}^{n+1} \| < \eta`, iteration is accepted, else compute a new iteration with the new time step.
4. Every :math:`N_\delta` iterations, or if :math:`\Delta t\not\in[\beta\Delta t^\star, \gamma\Delta t^\star]` (or if first iteration and :math:`\Delta t^\star` is not computed), compute new :math:`\Delta t^\star` from current :math:`\delta` (given by :cpp:func:`~ponio::splitting::strang::adaptive_strang::info` returned value which has a ``data`` member variable) and :math:`C_0` (given by :cpp:func:`ponio::splitting::strang::adaptive_strang::lipschitz_constant_estimate` method) with:

   .. math::

      \Delta t^\star \approx \frac{\delta C_\delta}{C_0},\qquad \delta C_\delta \Delta t^2 = \| \mathcal{S}^{\Delta t}(t^n, u^n) - \mathcal{S}_\delta^{\Delta t}(t^n, u^n) \|
5. If :math:`\Delta t \not\in [\beta\Delta t^\star, \gamma\Delta t^\star]`, assume :math:`\Delta t^\star=\Delta t` and compute a new :math:`\delta` with:

   .. math::

      \Delta t^\star \approx \frac{\delta C_\delta}{C_0}


----


IMEX methods with extended stability method
-------------------------------------------

The PIROCK method is introduce in :cite:`abdulle:2013`, the complete scheme is a IMEX scheme that allows for solving an equation of the following form (with 3 operators):

.. math::

  \dot{u} = F_R(u) + F_D(u) + F_A(u)

with :math:`F_R` a reaction operator (solved implicitly), :math:`F_D` a diffusion operator (solved explicitly by a modified ROCK2 method) and :math:`F_A` an advection operator (solved with an explicit RK3 method). The first implementation of PIROCK scheme in ponio makes to solve only reaction-diffusion problem (e.g. :math:`F_A = 0`), the second one to solve a complet reaction-diffusion-advection problem.

The method is divided into 5 steps:

1. Diffusion stabilization procedure (modified ROCK2 method)

   .. math::

      \begin{aligned}
        u^{(0)} &= u^n \\
        u^{(1)} &= u^n + \alpha \mu_1 \Delta t F_D(u^n) \\
        u^{(j)} &= \alpha \mu_j \Delta t F_D(u^{(j-1)}) - \nu_j u^{(j-1)} - \kappa_j u^{(j-2)},\quad j=2,\dots, s-2+\ell
      \end{aligned}
2. Finishing procedure for diffusion

   .. math::

      \begin{aligned}
         u^{*(s-1)} &= u^{(s-2)} + \sigma_\alpha \Delta t F_D(u^{(s-2)}) \\
         u^{*(s)} &= u^{*(s-1)} + \sigma_\alpha \Delta t F_D(u^{*(s-1)})
      \end{aligned}
3. Starting value for advection-reaction

   .. math::

      U = u^{(s-2+\ell)}
4. Finishing procedure for advection-reaction coupling

   .. math::

      \begin{aligned}
         u^{(s+1)} &= u^{(s-2+\ell)} + \gamma \Delta t F_R(u^{(s+1)}) \\
         u^{(s+2)} &= u^{(s-2+\ell)} + \beta \Delta t F_D(u^{(s+1)}) + \Delta t F_A(u^{(s+1)}) + (1-2\gamma)\Delta t F_R(u^{(s+1)}) + \gamma \Delta t F_R(u^{(s+2)}) \\
         u^{(s+3)} &= u^{(s-2+\ell)} + (1-2\gamma)\Delta t F_A(u^{(s+1)}) + (1-\gamma)\Delta t F_R(u^{(s+1)}) \\
         u^{(s+4)} &= u^{(s-2+\ell)} + \frac{1}{3}\Delta t F_A(u^{(s+1)}) \\
         u^{(s+5)} &= u^{(s-2+\ell)} + \frac{2\beta}{3} \Delta t F_D(u^{(s+1)}) + \frac{2}{3}\Delta t J_R^{-1} F_A(u^{(s+4)}) + \left(\frac{2}{3} - \gamma\right) \Delta t F_R(u^{(s+1)}) + \frac{2\gamma}{3}\Delta t F_R(u^{(s+2)})
      \end{aligned}
5. Computation of the integration step

   .. math::

      \begin{aligned}
         u^{n+1} =& u^{*(s)} \\
                     & - \sigma_\alpha\left( 1 - \frac{\tau_\alpha}{\sigma_\alpha^2}\right)\Delta t (F_D(u^{*(s-1)}) - F_D(u^{(s-2)})) \\
                     & + \frac{1}{4}\Delta t F_A(u^{(s+1)}) + \frac{3}{4}\Delta t F_A(u^{(s+5)}) \\
                     & + \frac{1}{2}\Delta t F_R(u^{(s+1)}) + \frac{1}{2}\Delta t F_R(u^{(s+2)}) \\
                     & + \frac{J_R^{-\ell}}{2-4\gamma}\Delta t (F_D(u^{(s+3)}) - F_D(u^{(s+1)}))
      \end{aligned}

   where :math:`\mu_j`, :math:`\nu_j`, :math:`\kappa_j` are the same coefficients as for the standard ROCK2 method and:

   * :math:`\gamma = 1-\frac{\sqrt{2}}{2}`
   * :math:`\beta = 1-2\alpha P'_{s-2+\ell}(0)` where :math:`P'_{s-2+\ell}` is the stability polynomial of the underlying ROCK2 method
   * :math:`J_R = I - \gamma \Delta t \frac{\partial F_R}{\partial u}(u^{(s-2+\ell)})`
   * :math:`\sigma_\alpha = \frac{1-\alpha}{2} + \alpha \sigma`, where :math:`\sigma` coming from the underlying ROCK2 method
   * :math:`\tau_\alpha = \frac{(\alpha-1)^2}{2} + 2\alpha(1-\alpha)\sigma + \alpha^2\tau`, where :math:`\sigma` and :math:`\tau` coming from the underlying ROCK2 method

For adaptive time step method, need to compute one error by operator and take the maximum value

.. math::

    \begin{aligned}
      err_D &= \sigma_\alpha \Delta t \left( 1 - \frac{\tau_\alpha}{\sigma_\alpha^2} \right) ( F_D( u^{*(s-1)} - u^{(s-2)} )) \\
      err_R &= \Delta t J_R^{-1}\left( \frac{1}{6}F_R(u^{(s+1)}) - \frac{1}{6}F_R(u^{(s+2)}) \right) \\
      err_A &= -\frac{3}{20}\Delta t F_A(u^{(s+1)}) + \frac{3}{10}\Delta t F_A(u^{(s+4)}) - \frac{3}{20}\Delta t F_A(u^{(s+5)})
    \end{aligned}

The error is computed with the maximum:

.. math::

  err = \max\left( \|err_D\|, \|err_R\|, \|err_A\|^{2/3} \right)

with the following error norm:

.. math::

  \|\cdot\| : err \mapsto \sum_i \left( \frac{err_i}{a_{tol} + \max(y^{n+1}_i, y^n_i) } \right)^2

where :math:`y^n` and :math:`y^{n+1}` are estimation of the solution respectively at time :math:`t^n` and after an iteration :math:`t^{n+1}=t^n+\Delta t`, and :math:`a_{tol}` and :math:`r_{tol}` respectively absolute and relative tolerances. This norm give a normalized error, so its compare to :math:`1` and new time step is computed as:

.. math::

  \Delta t_{new} = \sqrt{\frac{1}{err}} \Delta t

.. note::

  As in a lot of adaptive time step method, the new time step stays in :math:`\Delta t_{new} \in [0.5\Delta t, 2 \Delta t]` set and there is a safety factor, so the new time step is multiply by :math:`0.8`.

Still keep two free parameters :math:`\ell` and :math:`\alpha` given free to user, ponio provides two choices for this parameters:

* :math:`\ell=2` and :math:`\alpha=1`, in this case, if :math:`F_A=0` and :math:`F_R=0` we have the standard ROCK2 method;
* :math:`\ell=1` and :math:`\alpha = \frac{1}{2P'_{s-2+\ell}(0)}`, so :math:`\beta=0` to minimized computation cost.

.. seealso::

   See the :doc:`list of IMEX methods with extended stability method <../api/algorithm/list_alg_pirock>` in ponio.

----


Bibliography
------------

.. bibliography::
