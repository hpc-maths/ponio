List of IMEX stabilized methods
===============================


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


PIROCK for reaction-diffusion problem
-------------------------------------

In this section we present the implementation of PIROCK method for only reaction-diffusion problem (i.e. :math:`F_A = 0`).

.. doxygenclass:: ponio::runge_kutta::pirock::pirock_impl
   :project: ponio
   :members:

Helper functions
~~~~~~~~~~~~~~~~

User interface functions to build a PIROCK method.

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock(alpha_beta_computer_t&&, eig_computer_t&&, shampine_trick_caller_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock(alpha_beta_computer_t&&, eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock()
  :project: ponio


:math:`\ell=2` and :math:`\alpha = 1` case
""""""""""""""""""""""""""""""""""""""""""

Following functions are useful for to build a PIROCK method with :math:`\ell=2` and :math:`\alpha = 1` (with :cpp:class:`ponio::runge_kutta::pirock::alpha_fixed` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_a1(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_a1()
  :project: ponio


:math:`\ell=1` and :math:`\beta = 0` case
"""""""""""""""""""""""""""""""""""""""""

Following functions are useful for to build a PIROCK method with :math:`\ell=1` and :math:`\beta = 0` (with :cpp:class:`ponio::runge_kutta::pirock::beta_0` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_b0(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_b0()
  :project: ponio


PIROCK for reaction-diffusion-advection problem
-----------------------------------------------

In this section we present the implementation of complet PIROCK method (i.e. for reaction-diffusion-advection problem).

.. doxygenclass:: ponio::runge_kutta::pirock::pirock_RDA_impl
   :project: ponio
   :members:


Helper functions
~~~~~~~~~~~~~~~~

User interface functions to build a PIROCK method.

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA(alpha_beta_computer_t&&, eig_computer_t&&, shampine_trick_caller_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA(alpha_beta_computer_t&&, eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA()
  :project: ponio


:math:`\ell=2` and :math:`\alpha = 1` case
""""""""""""""""""""""""""""""""""""""""""

Following functions are useful for to build a PIROCK method with :math:`\ell=2` and :math:`\alpha = 1` (with :cpp:class:`ponio::runge_kutta::pirock::alpha_fixed` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_a1(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_a1()
  :project: ponio


:math:`\ell=1` and :math:`\beta = 0` case
"""""""""""""""""""""""""""""""""""""""""

Following functions are useful for to build a PIROCK method with :math:`\ell=1` and :math:`\beta = 0` (with :cpp:class:`ponio::runge_kutta::pirock::beta_0` computer).

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_b0(eig_computer_t&&)
  :project: ponio

.. doxygenfunction:: ponio::runge_kutta::pirock::pirock_RDA_b0()
  :project: ponio
