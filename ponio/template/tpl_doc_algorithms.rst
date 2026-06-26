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

Lawson methods was initialy propose into :cite:`lawson:1967`. Lawson methods are build to solve a problem with a linear and nonlinear part, to solve exactly the problem when the nonlinear part goes to zero. This class of problem can be write as

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

Some problems, like heat equation, require methods stabilized on the negative real axis. The ponio library provides a Runge-Kutta Chebyshev method of order 2, ROCK2 method (of order 2), ROCK4 method (of order 4) and a Runge-Kutta Legendre method of order 1 or 2.

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

.. doxygenfunction:: ponio::runge_kutta::rock::rock2(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::rock::rock2()
  :project: ponio


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

.. doxygenclass:: ponio::runge_kutta::rock::rock4_impl
   :project: ponio
   :members:

.. doxygenfunction:: ponio::runge_kutta::rock::rock4(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::rock::rock4()
  :project: ponio


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

.. doxygenclass:: ponio::runge_kutta::legendre::explicit_rkl1
   :project: ponio
   :members:


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

IMEX stabilized methods
-----------------------

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

.. warning::

  Adaptive time step method is not yet implemented for complet PIROCK :cpp:class:`ponio::runge_kutta::pirock::pirock_RDA_impl`.


Still keep two free parameters :math:`\ell` and :math:`\alpha` given free to user, ponio provides two choices for this parameters:

* :math:`\ell=2` and :math:`\alpha=1`, in this case, if :math:`F_A=0` and :math:`F_R=0` we have the standard ROCK2 method;
* :math:`\ell=1` and :math:`\alpha = \frac{1}{2P'_{s-2+\ell}(0)}`, so :math:`\beta=0` to minimized computation cost.

ponio provides two computers for :math:`\alpha` and :math:`\beta` values.

.. doxygenclass:: ponio::runge_kutta::pirock::alpha_fixed
   :project: ponio
   :members:

.. doxygenclass:: ponio::runge_kutta::pirock::beta_0
   :project: ponio
   :members:


PIROCK for reaction-diffusion problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we present the implementation of PIROCK method for only reaction-diffusion problem (i.e. :math:`F_A = 0`).

.. doxygenclass:: ponio::runge_kutta::pirock::pirock_impl
   :project: ponio
   :members:

Helper functions
""""""""""""""""

User interface functions to build a PIROCK method.

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock(alpha_beta_computer_t&&, eig_computer_t&&, shampine_trick_caller_t&&, value_t, value_t)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock(alpha_beta_computer_t&&, eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock()
  :project: ponio

Following functions are useful for to build a PIROCK method with :math:`\ell=2` and :math:`\alpha = 1` (with :cpp:class:`ponio::runge_kutta::pirock::alpha_fixed` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_a1(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_a1()
  :project: ponio

Following functions are useful for to build a PIROCK method with :math:`\ell=1` and :math:`\beta = 0` (with :cpp:class:`ponio::runge_kutta::pirock::beta_0` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_b0(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_b0()
  :project: ponio

With this functions you can initialize different kind of PIROCK method to solve your problem:

The mainly diffusive case with :math:`\alpha=1` and :math:`\ell=2` (to get a ROCK2 method when reaction becomes null) use :cpp:func:`ponio::runge_kutta::pirock::pirock_a1` or :cpp:func:`ponio::runge_kutta::pirock::pirock` with a :cpp:class:`ponio::runge_kutta::pirock::alpha_fixed` object, you can also add a `Shampine's trick` caller to get adaptive time step method.

.. code-block:: cpp

    auto pirock = ponio::runge_kutta::pirock::pirock_a1();

    // or

    auto adaptive_pirock = ponio::runge_kutta::pirock::pirock<2, true>(
      ponio::runge_kutta::pirock::alpha_fixed<double>( 1.0 ),
      eigmax_computer,                                        // [](auto& f, double tn, auto& un, double dt, auto& du_work)
      ponio::shampine_trick::shampine_trick<decltype( un )>() // Shampine's trick class
    );

The mainly reactive case with :math:`\beta=0` (or with less computational cost) you can also choose :math:`\ell = 1`, use :cpp:func:`ponio::runge_kutta::pirock::pirock_b0` or :cpp:func:`ponio::runge_kutta::pirock::pirock` with a :cpp:class:`ponio::runge_kutta::pirock::beta_0` object, you can also add a `Shampine's trick` caller to get adaptive time step method.

.. code-block:: cpp

    auto pirock = ponio::runge_kutta::pirock::pirock_b0();

    // or

    auto adaptive_pirock = ponio::runge_kutta::pirock::pirock<1, true>(
      ponio::runge_kutta::pirock::beta_0<double>(),
      eigmax_computer,                                        // [](auto& f, double tn, auto& un, double dt, auto& du_work)
      ponio::shampine_trick::shampine_trick<decltype( un )>() // Shampine's trick class
    );


PIROCK for reaction-diffusion-advection problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we present the implementation of complet PIROCK method (i.e. for reaction-diffusion-advection problem).

.. doxygenclass:: ponio::runge_kutta::pirock::pirock_RDA_impl
   :project: ponio
   :members:

Helper functions
""""""""""""""""

User interface functions to build a PIROCK method.

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA(alpha_beta_computer_t&&, eig_computer_t&&, shampine_trick_caller_t&&, value_t, value_t)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA(alpha_beta_computer_t&&, eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA()
  :project: ponio

Following functions are useful for to build a PIROCK method with :math:`\ell=2` and :math:`\alpha = 1` (with :cpp:class:`ponio::runge_kutta::pirock::alpha_fixed` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_a1(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_a1()
  :project: ponio

Following functions are useful for to build a PIROCK method with :math:`\ell=1` and :math:`\beta = 0` (with :cpp:class:`ponio::runge_kutta::pirock::beta_0` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_b0(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_b0()
  :project: ponio



Bibliography
------------

.. bibliography::
