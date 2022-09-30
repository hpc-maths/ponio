Algorithms
==========

List of all algorithms (skeleton of each numerical method to solve)

Runge-Kutta methods 
-------------------

Explicit methods
~~~~~~~~~~~~~~~~ 

.. doxygentypedef:: ode::butcher::euler
   :project: solver

.. doxygentypedef:: ode::butcher::explicit_euler_sub4
   :project: solver

.. doxygentypedef:: ode::butcher::rk6es
   :project: solver

.. doxygentypedef:: ode::butcher::rk_118
   :project: solver

.. doxygentypedef:: ode::butcher::rk_21
   :project: solver

.. doxygentypedef:: ode::butcher::rk_22_midpoint
   :project: solver

.. doxygentypedef:: ode::butcher::rk_22_ralston
   :project: solver

.. doxygentypedef:: ode::butcher::rk_32_best
   :project: solver

.. doxygentypedef:: ode::butcher::rk_33
   :project: solver

.. doxygentypedef:: ode::butcher::rk_33_233e
   :project: solver

.. doxygentypedef:: ode::butcher::rk_33_bogackishampine
   :project: solver

.. doxygentypedef:: ode::butcher::rk_33_heun
   :project: solver

.. doxygentypedef:: ode::butcher::rk_33_ralston
   :project: solver

.. doxygentypedef:: ode::butcher::rk_33_van_der_houwen
   :project: solver

.. doxygentypedef:: ode::butcher::rk_44
   :project: solver

.. doxygentypedef:: ode::butcher::rk_44_235j
   :project: solver

.. doxygentypedef:: ode::butcher::rk_44_38
   :project: solver

.. doxygentypedef:: ode::butcher::rk_44_ralston
   :project: solver

.. doxygentypedef:: ode::butcher::rk_65
   :project: solver

.. doxygentypedef:: ode::butcher::rk_65_236a
   :project: solver

.. doxygentypedef:: ode::butcher::rk_76
   :project: solver

.. doxygentypedef:: ode::butcher::rk_86
   :project: solver

.. doxygentypedef:: ode::butcher::rk_nssp_21
   :project: solver

.. doxygentypedef:: ode::butcher::rk_nssp_32
   :project: solver

.. doxygentypedef:: ode::butcher::rk_nssp_33
   :project: solver

.. doxygentypedef:: ode::butcher::rk_nssp_53
   :project: solver

.. doxygentypedef:: ode::butcher::rk_spp_43
   :project: solver

.. doxygentypedef:: ode::butcher::rk_ssp_22_heun
   :project: solver

.. doxygentypedef:: ode::butcher::rk_ssp_32
   :project: solver

.. doxygentypedef:: ode::butcher::rk_ssp_33
   :project: solver

.. doxygentypedef:: ode::butcher::rk_ssp_42
   :project: solver

.. doxygentypedef:: ode::butcher::rk_ssp_53
   :project: solver

.. doxygentypedef:: ode::butcher::rk_ssp_54
   :project: solver

Embedded methods
~~~~~~~~~~~~~~~~

.. doxygentypedef:: ode::butcher::rk54_6m
   :project: solver

.. doxygentypedef:: ode::butcher::rk54_7m
   :project: solver

.. doxygentypedef:: ode::butcher::rk54_7s
   :project: solver


Lawson methods 
--------------

This Lawson methods have a underlying Runge-Kutta method, they have the same name just prefixed by ``l``.

Lawson methods are build to solve a problem like:

.. math::

  \dot{v} = Lu + N(t,v)

First, we introduce the change of variable

.. math::

   u(t) = e^{-Lt}v(t)

which yields the equation

.. math::

   \dot{u}(t) = -Le^{-Lt}v(t) + e^{-Lt}\dot{v}(t)

which can be rewrite in term of :math:`u` as

.. math::

   \dot{u}(t) = e^{-Lt}N(t,e^{Lt}u) = \tilde{N}(t,u)

We solve this equation with a classical Runge-Kutta method RK(:math:`s`,:math:`n`) with :math:`s` stages and of order :math:`n`. We rewrite the scheme in term of the variable :math:`v`, which yields the following scheme

.. math::

   \begin{aligned}
      k_i &= e^{-c_i\Delta tL}N(t^n+c_i\Delta t, e^{c_i\Delta tL}\left( v^n+\sum_{j}a_{i,j}k_j \right))\\
      v^{n+1} &= e^{\Delta tL}\left( v^n + \Delta t \sum_i b_i k_i \right)
   \end{aligned}

Explicit methods
~~~~~~~~~~~~~~~~ 

.. doxygenfunction:: ode::butcher::leuler
   :project: solver

.. doxygenfunction:: ode::butcher::lexplicit_euler_sub4
   :project: solver

.. doxygenfunction:: ode::butcher::lrk6es
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_118
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_21
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_22_midpoint
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_22_ralston
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_32_best
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_33
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_33_233e
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_33_bogackishampine
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_33_heun
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_33_ralston
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_33_van_der_houwen
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_44
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_44_235j
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_44_38
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_44_ralston
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_65
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_65_236a
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_76
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_86
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_nssp_21
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_nssp_32
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_nssp_33
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_nssp_53
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_spp_43
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_ssp_22_heun
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_ssp_32
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_ssp_33
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_ssp_42
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_ssp_53
   :project: solver

.. doxygenfunction:: ode::butcher::lrk_ssp_54
   :project: solver


Embedded methods
~~~~~~~~~~~~~~~~

.. doxygenfunction:: ode::butcher::lrk54_6m
   :project: solver

.. doxygenfunction:: ode::butcher::lrk54_7m
   :project: solver

.. doxygenfunction:: ode::butcher::lrk54_7s
   :project: solver

Runge-Kutta Chebyshev methods
-----------------------------

.. doxygenclass:: ode::butcher::chebyshev::explicit_rkc2
   :project: solver
   :members:

