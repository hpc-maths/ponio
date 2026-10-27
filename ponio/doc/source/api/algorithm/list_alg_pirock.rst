List of IMEX stabilized methods
===============================

The PIROCK method is introduced in :cite:`abdulle:2013`. It is a second-order partitioned IMEX scheme for problems of the form

.. math::

   \dot{u} = F_R(u) + F_D(u) + F_A(u),

where :math:`F_R` is a reaction operator treated implicitly, :math:`F_D` is a diffusion operator stabilized with a modified ROCK2 method and :math:`F_A` is an advection operator treated explicitly.

The ponio implementation applies PIROCK directly to the complete state. This differs from the reference Fortran implementation, where the reaction solves are performed locally by blocks of size ``npdes``. Three variants are available in ponio: reaction-diffusion (RD, :math:`F_A=0`), diffusion-advection (DA, :math:`F_R=0`) and the complete reaction-diffusion-advection (RDA) problem.

Common diffusion stages
-----------------------

All three variants use the same stabilized diffusion procedure. With :math:`s`
the ROCK2 number of stages,

.. math::

   \begin{aligned}
      u^{(0)} &= u^n, \\
      u^{(1)} &= u^n + \alpha\mu_1\Delta t F_D(u^n), \\
      u^{(j)} &= \alpha\mu_j\Delta t F_D(u^{(j-1)})
               - \nu_j u^{(j-1)} - \kappa_j u^{(j-2)}.
   \end{aligned}

The recurrence is first performed for

.. math::

   j=2,\dots,s-2,

using the standard ROCK2 recurrence coefficients. The state
:math:`u^{(s-2)}` is then kept for the diffusion finishing procedure.

The recurrence is subsequently extended by :math:`\ell` supplementary stages,

.. math::

   j=s-1,\dots,s-2+\ell,

using the corresponding supplementary ROCK2 coefficients. These stages produce

.. math::

   U = u^{(s-2+\ell)},

which is the starting value for the partitioned finishing procedure. Thus,
:math:`u^{(s-2)}` and :math:`U` play distinct roles: the former is used to
construct the diffusion finishing stages, while the latter includes the
:math:`\ell` supplementary stages and is used for the subsequent coupling with
the other partitions.

Starting from :math:`u^{(s-2)}`, the first diffusion finishing stage is

.. math::

   u^{*(s-1)}
   =
   u^{(s-2)}
   +
   \sigma_\alpha\Delta t F_D(u^{(s-2)}).

A second diffusion evaluation at :math:`u^{*(s-1)}` provides the embedded
diffusion error estimator

.. math::

   err_D
   =
   \sigma_\alpha
   \left(
      1-\frac{\tau_\alpha}{\sigma_\alpha^2}
   \right)
   \Delta t
   \left(
      F_D(u^{*(s-1)})-F_D(u^{(s-2)})
   \right).

The second finishing stage then gives the completed diffusion contribution

.. math::

   u_D
   =
   u^{*(s-1)}
   +
   \sigma_\alpha\Delta t F_D(u^{*(s-1)})
   -
   err_D.

Hence the stabilized diffusion construction produces two quantities used later
in PIROCK: :math:`U=u^{(s-2+\ell)}`, which enters the partitioned finishing
procedure, and :math:`u_D`, which contains the completed diffusion contribution.

The ROCK2 degree is selected dynamically at each step from the diffusion
stability requirement. Ponio uses the ROCK2 stability constants :math:`0.432`
for :math:`\ell=1` and :math:`0.811` for :math:`\ell=2`, and sets
:math:`s=mdeg+2` from a supported ROCK2 degree.

The coefficients are defined by

.. math::

   \gamma = 1-\frac{\sqrt{2}}{2},
   \qquad
   \beta = 1-2\alpha P'_{s-2+\ell}(0),

and

.. math::

   \sigma_\alpha = \frac{1-\alpha}{2}+\alpha\sigma,
   \qquad
   \tau_\alpha
   =
   \frac{(\alpha-1)^2}{2}
   +2\alpha(1-\alpha)\sigma
   +\alpha^2\tau,

where :math:`\mu_j`, :math:`\nu_j`, :math:`\kappa_j`, :math:`\sigma` and
:math:`\tau` are the ROCK2 coefficients.

Reaction-diffusion finishing
----------------------------

For :math:`F_A=0`, only the reaction-diffusion coupling is required:

.. math::

   \begin{aligned}
      u^{(s+1)} &= U + \gamma\Delta t F_R(u^{(s+1)}), \\
      u^{(s+2)} &= U + \beta\Delta t F_D(u^{(s+1)})
                    + (1-2\gamma)\Delta t F_R(u^{(s+1)})
                    + \gamma\Delta t F_R(u^{(s+2)}), \\
      u^{(s+3)} &= U + (1-\gamma)\Delta t F_R(u^{(s+1)}).
   \end{aligned}

With

.. math::

   J_R = I - \gamma\Delta t\frac{\partial F_R}{\partial u}(U),

one step is completed by

.. math::

   u^{n+1} = u_D
      + \frac{\Delta t}{2}\left(F_R(u^{(s+1)})+F_R(u^{(s+2)})\right)
      + \frac{\Delta t}{2-4\gamma}J_R^{-\ell}
        \left(F_D(u^{(s+3)})-F_D(u^{(s+1)})\right).

For Jacobian-based reaction problems, the implicit stages use a simplified Newton iteration following the ``ieuler`` routine of the reference Fortran code (lines 990--1054). The same reaction Jacobian and its factorization are reused for the first five Newton corrections. If convergence is not reached, the Jacobian is refreshed before correction 6 and before each following correction. In ponio this solve acts on the complete state, instead of the local ``npdes`` blocks used by the Fortran implementation.

The factor :math:`J_R^{-\ell}` is the Shampine damping used in the diffusion-reaction coupling. The same idea is used for the reaction error estimator with :math:`J_R^{-1}`. When the reaction matrix is already factorized by the implicit stages, ponio reuses this factorization for the Shampine solves. The corresponding reaction-error damping can be seen in the Fortran routines around ``dampreaction`` (lines 945--986).

Diffusion-advection finishing
-----------------------------

For :math:`F_R=0`, no implicit solve and no Shampine damping are needed. Setting :math:`K=U`, the remaining stages are

.. math::

   \begin{aligned}
      u^{(s+3)} &= K + (1-2\gamma)\Delta t F_A(K), \\
      u^{(s+4)} &= K + \frac{1}{3}\Delta t F_A(K), \\
      u^{(s+5)} &= K + \frac{2\beta}{3}\Delta t F_D(K)
                    + \frac{2}{3}\Delta t F_A(u^{(s+4)}),
   \end{aligned}

and

.. math::

   \begin{aligned}
      u^{n+1} = u_D
        &+ \frac{1}{4}\Delta t F_A(K)
         + \frac{3}{4}\Delta t F_A(u^{(s+5)}) \\
        &+ \frac{\Delta t}{2-4\gamma}
           \left(F_D(u^{(s+3)})-F_D(K)\right).
   \end{aligned}

Reaction-diffusion-advection finishing
--------------------------------------

For the complete RDA problem, the finishing stages are

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

The final update is

.. math::

   \begin{aligned}
      u^{n+1} = u_D
        &+ \frac{1}{4}\Delta t F_A(u^{(s+1)})
         + \frac{3}{4}\Delta t F_A(u^{(s+5)}) \\
        &+ \frac{1}{2}\Delta t F_R(u^{(s+1)})
         + \frac{1}{2}\Delta t F_R(u^{(s+2)}) \\
        &+ \frac{\Delta t}{2-4\gamma}J_R^{-\ell}
           \left(F_D(u^{(s+3)})-F_D(u^{(s+1)})\right).
   \end{aligned}

The implicit reaction stages and the Shampine solves use the same strategy as in the RD variant.

Adaptive time step
------------------

The embedded PIROCK variants compute one estimator for each active operator:

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

Ponio stores the normalized squared RMS error

.. math::

   E(e) = \frac{1}{N}\sum_{i=1}^{N}
      \left(
         \frac{|e_i|}
              {a_{tol}+r_{tol}\max(|u_i^n|,|u_i^{n+1}|)}
      \right)^2.

The error used by the controller is therefore

.. math::

   \begin{aligned}
      E_{RD}  &= \max(E(err_D),E(err_R)), \\
      E_{DA}  &= \max(E(err_D),E(err_A)^{2/3}), \\
      E_{RDA} &= \max(E(err_D),E(err_R),E(err_A)^{2/3}).
   \end{aligned}

A step is accepted when :math:`E<1`. Since the stored quantity is an RMS error squared, the raw step-size factor is :math:`E^{-1/4}`. The controller then applies the safety factor, growth bounds and, when diffusion dominates, the memory correction. A rejected attempt keeps :math:`(t^n,u^n)` unchanged and retries with a smaller step; the diffusion evaluation at :math:`(t^n,u^n)` is reused on the retry. The diffusion estimator can also reject a step before the more expensive reaction/advection finishing stages are computed.

This follows the error-control strategy of the reference Fortran code, in particular the controller in lines 520--615 and the normalized error estimates around lines 945--970.

The two standard parameter choices provided by ponio are

* :math:`\ell=2` and :math:`\alpha=1`, with :cpp:class:`ponio::runge_kutta::pirock::alpha_fixed`;
* :math:`\ell=1` and :math:`\alpha = \frac{1}{2P'_{s-2+\ell}(0)}`, hence :math:`\beta=0`, with :cpp:class:`ponio::runge_kutta::pirock::beta_0`.

PIROCK for reaction-diffusion problem
-------------------------------------

+ **name:** PIROCK for reaction-diffusion problem
+ **label in ponio:** :cpp:class:`ponio::runge_kutta::pirock::pirock_impl`
+ **helper function in ponio:** :cpp:func:`ponio::runge_kutta::pirock::pirock`
+ **stages:** dynamic number of stages
+ **order:** 2

Following the choice of parameters, ponio provides two fixed-step helper functions:

+ :math:`\ell=2` and :math:`\alpha=1`, with :cpp:func:`ponio::runge_kutta::pirock::pirock_a1`
+ :math:`\ell=1` and :math:`\beta=0`, with :cpp:func:`ponio::runge_kutta::pirock::pirock_b0`

PIROCK for diffusion-advection problem
--------------------------------------

+ **name:** PIROCK for diffusion-advection problem
+ **label in ponio:** :cpp:class:`ponio::runge_kutta::pirock::pirock_DA_impl`
+ **helper function in ponio:** :cpp:func:`ponio::runge_kutta::pirock::pirock_DA`
+ **stages:** dynamic number of stages
+ **order:** 2

Following the choice of parameters, ponio provides two fixed-step helper functions:

+ :math:`\ell=2` and :math:`\alpha=1`, with :cpp:func:`ponio::runge_kutta::pirock::pirock_DA_a1`
+ :math:`\ell=1` and :math:`\beta=0`, with :cpp:func:`ponio::runge_kutta::pirock::pirock_DA_b0`

PIROCK for reaction-diffusion-advection problem
-----------------------------------------------

+ **name:** PIROCK for reaction-diffusion-advection problem
+ **label in ponio:** :cpp:class:`ponio::runge_kutta::pirock::pirock_RDA_impl`
+ **helper function in ponio:** :cpp:func:`ponio::runge_kutta::pirock::pirock_RDA`
+ **stages:** dynamic number of stages
+ **order:** 2

Following the choice of parameters, ponio provides two fixed-step helper functions:

+ :math:`\ell=2` and :math:`\alpha=1`, with :cpp:func:`ponio::runge_kutta::pirock::pirock_RDA_a1`
+ :math:`\ell=1` and :math:`\beta=0`, with :cpp:func:`ponio::runge_kutta::pirock::pirock_RDA_b0`
