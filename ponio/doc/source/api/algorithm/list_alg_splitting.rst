List of splitting methods
=========================

Some methods to solve a problem of the form:

.. math::

  \dot{u} = \sum_i f_i(t,u)

with the initial condition :math:`u(t=0)=u_0`. We note with :math:`\phi_{\tau}^{[f_i]}(t^n,\tilde{u}^n)` the solution at time :math:`t^n+\tau` of the subproblem :math:`\dot{u}=f_i(t,u)` with the initial condition :math:`u(t^n)=\tilde{u}^n`.

Lie splitting method
--------------------

In Lie splitting method, the solution is computed as:

.. math::

   u^{n+1} = \mathcal{L}^{\Delta t}(t^n, u^n) = \phi_{\Delta t}^{[f_1]}\circ \cdots \circ \phi_{\Delta t}^{[f_k]} (t^n,u^n)

+ **name:** Lie splitting method
+ **label in ponio:** :cpp:class:`ponio::splitting::lie::lie`
+ **helper function in ponio:** :cpp:class:`ponio::splitting::lie::make_lie_tuple`
+ **order:** 1


Strang splitting method
-----------------------

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \mathcal{S}^{\Delta t}(t^n, u^n) = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}
              \circ \phi_{\Delta t}^{[f_k]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{k-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

+ **name:** Strang splitting method
+ **label in ponio:** :cpp:class:`ponio::splitting::strang::strang`
+ **helper function in ponio:** :cpp:class:`ponio::splitting::strang::make_strang_tuple`
+ **order:** 2


Adaptive time step Strang splitting method
------------------------------------------

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
4. Every :math:`N_\delta` iterations, or if :math:`\Delta t\not\in[\beta\Delta t^\star, \gamma\Delta t^\star]` (or if first iteration and :math:`\Delta t^\star` is not computed), compute new :math:`\Delta t^\star` from current :math:`\delta` (given by :cpp:func:`~ponio::splitting::strang::adaptive_strang::info` returned value which has a ``data`` member variable) and :math:`C_0` (given by :cpp:func:`~ponio::splitting::strang::adaptive_strang::lipschitz_constant_estimate` method) with:

   .. math::

      \Delta t^\star \approx \frac{\delta C_\delta}{C_0},\qquad \delta C_\delta \Delta t^2 = \| \mathcal{S}^{\Delta t}(t^n, u^n) - \mathcal{S}_\delta^{\Delta t}(t^n, u^n) \|
5. If :math:`\Delta t \not\in [\beta\Delta t^\star, \gamma\Delta t^\star]`, assume :math:`\Delta t^\star=\Delta t` and compute a new :math:`\delta` with:

   .. math::

      \Delta t^\star \approx \frac{\delta C_\delta}{C_0}


+ **name:** adaptive Strang splitting method
+ **label in ponio:** :cpp:class:`ponio::splitting::strang::adaptive_strang`
+ **helper function in ponio:** :cpp:class:`ponio::splitting::strang::make_adaptive_strang_tuple`
+ **order:** 2
