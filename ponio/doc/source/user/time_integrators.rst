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

We will solve this problem with different kind of methods just to present how to use them.

.. raw:: html
    :file: ../_static/cpp/ch_all.html


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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/erk.cxx
  :language: cpp
  :lines: 1-4

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/erk.cxx
  :language: cpp
  :lines: 6

Embedded methods
~~~~~~~~~~~~~~~~

Some Butcher tableaus provide an extra :math:`b` line:

.. math::

   \begin{array}{c|c}
      c & A \\
      \hline
        & b^\top \\
      \hline
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

Most common embedded Runge-Kutta methods come from :footcite:t:`dormand:1980` and :footcite:t:`prince:1981`, that why tey are sometime call Dormand-Prince methods.

.. seealso::

   See the :doc:`list of embedded Runge-Kutta methods <../api/algorithm/list_alg_dp>` in ponio.

The most common embedded methods are explicit methods, so we can only define a the function :math:`f` as a lambda as

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/dp.cxx
  :language: cpp
  :lines: 1-4

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/dp.cxx
  :language: cpp
  :lines: 6

Diagonal implicit methods
~~~~~~~~~~~~~~~~~~~~~~~~~

When the matrix :math:`A` is lower triangular with a diagonal, the Runge-Kutta method is called diagonal implicit (or DIRK). You have to provide a Jacobian function that returns the Jacobian matrix in point :math:`(t, u)` (see :cpp:class:`ponio::implicit_problem`). You can also provide an operator base definition (see :cpp:class:`ponio::implicit_operator_problem`).

.. seealso::

   See the :doc:`list of diagonal implicit Runge-Kutta methods <../api/algorithm/list_alg_dirk>` in ponio.

For an diagonal implicit Runge-Kutta method, we need the function :math:`f` and also its Jacobian function :math:`\mathrm{d}f`, and store them into a :cpp:class:`ponio::implicit_problem`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/dirk.cxx
  :language: cpp
  :lines: 1-10

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/dirk.cxx
  :language: cpp
  :lines: 12


Lawson methods
~~~~~~~~~~~~~~

Lawson methods was initially propose into :footcite:t:`lawson:1967`. Lawson methods are build to solve a problem with a linear and nonlinear part, to solve exactly the problem when the nonlinear part goes to zero. This class of problem can be write as

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

For a problem split into a linear and nonlinear part, we need to define linear part as a scalar or a matrix and the nonlinear part as a function, and store them into a :cpp:class:`ponio::lawson_problem`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/lrk.cxx
  :language: cpp
  :lines: 1-7

Because sometime we want to define exponential function with a Pade approximant, we need to specify an exponential function.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/lrk.cxx
  :language: cpp
  :lines: 9-12

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/lrk.cxx
  :language: cpp
  :lines: 14

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

For a problem split into a linear and nonlinear part, we need to define linear part as a scalar or a matrix and the nonlinear part as a function, and store them into a :cpp:class:`ponio::lawson_problem`.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/exprk.cxx
  :language: cpp
  :lines: 1-7

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/exprk.cxx
  :language: cpp
  :lines: 9


Additive Runge-Kutta methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An additive Runge-Kutta method is defined in ponio by its pair of Butcher tableaus

.. math::

   \begin{array}{c|c}
      c & A \\
      \hline
        & b^\top
   \end{array}
   \quad
   \begin{array}{c|c}
      \tilde{c} & \tilde{A} \\
      \hline
        & \tilde{b}^\top
   \end{array}

which is stored as a JSON file in ``database`` folder, and where :math:`(A, b, c)` is an explicit Runge-Kutta method and :math:`(\tilde{A}, \tilde{b}, \tilde{c})` is a diagonal implicit Runge-Kutta method.

Additive Runge-Kutta methods are useful to solve a problem like

.. math::

   \dot{u} = f_e(t, u) + f_i(t,u)

The method, with the previous pair of Butcher tableaus, reads as

.. math::

   \begin{aligned}
      u^{(i)} &= u^n + \Delta t \sum_j a_{ij} k_j + \Delta t \sum_\ell \tilde{a}_{i\ell} \tilde{k}_\ell, \quad i = 1, \dots, s \\
      k_i         &= f_e(t^n + c_i\Delta t, u^{(i)}) \\
      \tilde{k}_i &= f_i(t^n + \tilde{c}_i\Delta t, u^{(i)}) \\
      u^{n+1} &= u^n + \Delta t \sum_i b_i k_i  + \Delta t \sum_j \tilde{b}_j k_j
   \end{aligned}

.. seealso::

   See the :doc:`list of additive Runge-Kutta methods <../api/algorithm/list_alg_ark>` in ponio.

For a problem split into a part which will be solved explicitly and a part which will be solved implicitly (you should provide the Jacobian matrix), and store them into a :cpp:class:`ponio::imex_problem`:

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/ark.cxx
  :language: cpp
  :lines: 1-14

And next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/ark.cxx
  :language: cpp
  :lines: 16

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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rkc.cxx
  :language: cpp
  :lines: 1-4

We need to specify how many stages we want (we choose 5), and next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rkc.cxx
  :language: cpp
  :lines: 6


ROCK2 method
~~~~~~~~~~~~

We write the method ROCK2 presented in :footcite:t:`abdulle:2001`. The algorithm of ROCK2 method is the following:

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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rock2.cxx
  :language: cpp
  :lines: 1-4

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rock2.cxx
  :language: cpp
  :lines: 6

.. hint::

   The ROCK2 method computes dynamically the number of stages with power method to estimate the spectral radius of your operator, but you can provide an other estimator of spectral radius if needed.

ROCK4 method
~~~~~~~~~~~~

We write the method ROCK2 presented in :footcite:t:`abdulle:2002`. The algorithm of ROCK4 method is the following:

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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rock4.cxx
  :language: cpp
  :lines: 1-4

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rock4.cxx
  :language: cpp
  :lines: 6

.. hint::

   The ROCK4 method computes dynamically the number of stages with power method to estimate the spectral radius of your operator, but you can provide an other estimator of spectral radius if needed.


Runge-Kutta Legendre method
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An other way to get a stabilized Runge-Kutta method is to use Legendre polynomials, we follow presentation in :footcite:t:`meyer:2014`.

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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rkl1.cxx
  :language: cpp
  :lines: 1-4

We need to specify how many stages we want (we choose 5), and next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rkl1.cxx
  :language: cpp
  :lines: 6

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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rkl2.cxx
  :language: cpp
  :lines: 1-4

We need to specify how many stages we want (we choose 5), and next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/rkl2.cxx
  :language: cpp
  :lines: 6

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

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/lie.cxx
  :language: cpp
  :lines: 1-10

We also define the Lie splitting method and the method to solve each substep. In this example we use two RK(3,3) methods, and we specify the time step for subcycle

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/lie.cxx
  :language: cpp
  :lines: 12-13

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/lie.cxx
  :language: cpp
  :lines: 15


Strang Splitting method
~~~~~~~~~~~~~~~~~~~~~~~

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \mathcal{S}^{\Delta t}(t^n, u^n) = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}
              \circ \phi_{\Delta t}^{[f_k]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

For splitting methods in ponio, we first need to define all sub-problem (here only :math:`f_1` and :math:`f_2`) and store them into a :cpp:class:`ponio::problem`

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/strang.cxx
  :language: cpp
  :lines: 1-10

We also define the Strang splitting method and the method to solve each substep. In this example we use two RK(3,3) methods, and we specify the time step for subcycle

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/strang.cxx
  :language: cpp
  :lines: 12-13

Next call the :cpp:func:`ponio::solve` function with

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/strang.cxx
  :language: cpp
  :lines: 15

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

The difficulty in adaptive time step Strang splitting method is to find a good shift coefficient :math:`\delta` for a given tolerance :math:`\eta`. The complete strategy of an iteration provides in :footcite:t:`descombes:2015` is

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

The PIROCK method is introduced in :footcite:t:`abdulle:2013`. It is a second-order partitioned IMEX scheme for

.. math::

   \dot{u} = F_R(u) + F_D(u) + F_A(u),

where :math:`F_R` is treated implicitly, :math:`F_D` is stabilized by a modified ROCK2 method and :math:`F_A` is treated explicitly. Ponio provides three variants: reaction-diffusion (RD, :math:`F_A=0`), diffusion-advection (DA, :math:`F_R=0`) and reaction-diffusion-advection (RDA).

The complete method can be written in five steps:

1. Diffusion stabilization procedure (modified ROCK2 method)

   .. math::

      \begin{aligned}
         u^{(0)} &= u^n, \\
         u^{(1)} &= u^n + \alpha\mu_1\Delta t F_D(u^n), \\
         u^{(j)} &= \alpha\mu_j\Delta t F_D(u^{(j-1)})
                    - \nu_j u^{(j-1)} - \kappa_j u^{(j-2)},
                    \quad j=2,\dots,s-2+\ell.
      \end{aligned}

2. Finishing procedure for diffusion

   .. math::

      \begin{aligned}
         u^{*(s-1)} &= u^{(s-2)} + \sigma_\alpha\Delta t F_D(u^{(s-2)}), \\
         u^{*(s)} &= u^{*(s-1)} + \sigma_\alpha\Delta t F_D(u^{*(s-1)}).
      \end{aligned}

3. Starting value for advection-reaction

   .. math::

      U = u^{(s-2+\ell)}.

4. Finishing procedure for advection-reaction coupling

   .. math::

      \begin{aligned}
         u^{(s+1)} &= U + \gamma\Delta t F_R(u^{(s+1)}), \\
         u^{(s+2)} &= U + \beta\Delta t F_D(u^{(s+1)}) + \Delta t F_A(u^{(s+1)})
                       + (1-2\gamma)\Delta t F_R(u^{(s+1)})
                       + \gamma\Delta t F_R(u^{(s+2)}), \\
         u^{(s+3)} &= U + (1-2\gamma)\Delta t F_A(u^{(s+1)})
                       + (1-\gamma)\Delta t F_R(u^{(s+1)}), \\
         u^{(s+4)} &= U + \frac{1}{3}\Delta t F_A(u^{(s+1)}), \\
         u^{(s+5)} &= U + \frac{2\beta}{3}\Delta t F_D(u^{(s+1)})
                       + \frac{2}{3}\Delta t J_R^{-1}F_A(u^{(s+4)})
                       + \left(\frac{2}{3}-\gamma\right)\Delta t F_R(u^{(s+1)})
                       + \frac{2\gamma}{3}\Delta t F_R(u^{(s+2)}).
      \end{aligned}

5. Computation of the integration step

   .. math::

      \begin{aligned}
         u^{n+1} =& u^{*(s)} \\
                    &- \sigma_\alpha\left(1-\frac{\tau_\alpha}{\sigma_\alpha^2}\right)\Delta t
                       \left(F_D(u^{*(s-1)})-F_D(u^{(s-2)})\right) \\
                    &+ \frac{1}{4}\Delta t F_A(u^{(s+1)})
                     + \frac{3}{4}\Delta t F_A(u^{(s+5)}) \\
                    &+ \frac{1}{2}\Delta t F_R(u^{(s+1)})
                     + \frac{1}{2}\Delta t F_R(u^{(s+2)}) \\
                    &+ \frac{\Delta t}{2-4\gamma}J_R^{-\ell}
                       \left(F_D(u^{(s+3)})-F_D(u^{(s+1)})\right).
      \end{aligned}

where

* :math:`\gamma = 1-\frac{\sqrt{2}}{2}`;
* :math:`\beta = 1-2\alpha P'_{s-2+\ell}(0)`;
* :math:`J_R = I-\gamma\Delta t\frac{\partial F_R}{\partial u}(U)`;
* :math:`\sigma_\alpha = \frac{1-\alpha}{2}+\alpha\sigma`;
* :math:`\tau_\alpha = \frac{(\alpha-1)^2}{2}+2\alpha(1-\alpha)\sigma+\alpha^2\tau`.

The RD and DA variants are obtained directly by setting :math:`F_A=0` and :math:`F_R=0`, respectively. In particular, DA does not require an implicit reaction solve or Shampine damping.

For adaptive time stepping, one estimator is computed for each active operator:

.. math::

   \begin{aligned}
      err_D &= \sigma_\alpha\left(1-\frac{\tau_\alpha}{\sigma_\alpha^2}\right)\Delta t
               \left(F_D(u^{*(s-1)})-F_D(u^{(s-2)})\right), \\
      err_R &= \frac{\Delta t}{6}J_R^{-1}
               \left(F_R(u^{(s+2)})-F_R(u^{(s+1)})\right), \\
      err_A &= -\frac{3}{20}\Delta t F_A(u^{(s+1)})
                + \frac{3}{10}\Delta t F_A(u^{(s+4)})
                - \frac{3}{20}\Delta t F_A(u^{(s+5)}).
   \end{aligned}

Ponio uses the normalized squared RMS error

.. math::

   E(e) = \frac{1}{N}\sum_{i=1}^{N}
      \left(
         \frac{|e_i|}
              {a_{tol}+r_{tol}\max(|u_i^n|,|u_i^{n+1}|)}
      \right)^2.

The controller uses :math:`\max(E(err_D),E(err_R),E(err_A)^{2/3})` for RDA, with the inactive terms removed for RD and DA. A step is accepted when this value is smaller than one. The step-size controller also applies a safety factor, growth limits and a memory correction when the diffusion error dominates.

Two standard parameter choices are provided:

* :math:`\ell=2` and :math:`\alpha=1`;
* :math:`\ell=1` and :math:`\alpha=\frac{1}{2P'_{s-2+\ell}(0)}`, hence :math:`\beta=0`.

.. seealso::

   See the :doc:`list of IMEX methods with extended stability method <../api/algorithm/list_alg_pirock>` in ponio for the detailed RD, DA and RDA algorithms.

Using PIROCK in ponio
~~~~~~~~~~~~~~~~~~~~~

For a reaction-diffusion problem, define :math:`F_D`, :math:`F_R` and the reaction Jacobian :math:`\partial_u F_R` with :cpp:func:`ponio::make_imex_jacobian_problem`:

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/pirock.cxx
   :language: cpp
   :lines: 1-14

The fixed-step helpers are :cpp:func:`ponio::runge_kutta::pirock::pirock_a1` for :math:`\ell=2`, :math:`\alpha=1`, and :cpp:func:`ponio::runge_kutta::pirock::pirock_b0` for :math:`\ell=1`, :math:`\beta=0`. For example,

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_all/pirock.cxx
   :language: cpp
   :lines: 16

For a diffusion-advection problem, the problem contains the diffusion operator first and the advection operator second:

.. code-block:: cpp

   auto pb = ponio::make_problem( fd, fa );
   auto method = ponio::runge_kutta::pirock::pirock_DA_b0();

   ponio::solve( pb, method, y_0, t_span, dt, observer );

The alternative fixed-step helper is :cpp:func:`ponio::runge_kutta::pirock::pirock_DA_a1`.

For a reaction-diffusion-advection problem, the reaction problem is stored first, followed by diffusion and advection:

.. code-block:: cpp

   auto reaction = ponio::make_implicit_problem( fr, dfr );
   auto pb = ponio::make_problem( reaction, fd, fa );
   auto method = ponio::runge_kutta::pirock::pirock_RDA_b0();

   ponio::solve( pb, method, y_0, t_span, dt, observer );

The alternative fixed-step helper is :cpp:func:`ponio::runge_kutta::pirock::pirock_RDA_a1`. Operator-based reaction problems can be built in the same way with the corresponding ponio operator-problem factories.

Adaptive time step and Shampine's trick
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The adaptive variants are selected with the ``is_embedded`` template parameter and tolerances are set with ``abs_tol`` and ``rel_tol``. For RD and RDA, the adaptive reaction estimator uses Shampine's trick; the Shampine caller is therefore passed to the method:

.. code-block:: cpp

   auto method = ponio::runge_kutta::pirock::pirock<1, true>(
       ponio::runge_kutta::pirock::beta_0<double>(),
       eig_computer,
       ponio::shampine_trick::shampine_trick<decltype( y_0 )>() );

   method.abs_tol( 1e-8 ).rel_tol( 1e-6 );

For RDA, use :cpp:func:`ponio::runge_kutta::pirock::pirock_RDA` in the same way. Shampine's trick can also be enabled in fixed-step RD/RDA by passing the same caller with ``is_embedded=false``.

DA has no reaction operator and therefore no Shampine step. Its adaptive variant is built directly with :cpp:func:`ponio::runge_kutta::pirock::pirock_DA`:

.. code-block:: cpp

   auto method = ponio::runge_kutta::pirock::pirock_DA<1, true>(
       ponio::runge_kutta::pirock::beta_0<double>(), eig_computer );

   method.abs_tol( 1e-8 ).rel_tol( 1e-6 );


----

Bibliography
------------

.. footbibliography::
