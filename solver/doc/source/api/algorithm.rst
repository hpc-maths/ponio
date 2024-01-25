Algorithms
==========

List of all algorithms (skeleton of each numerical method to solve)

Runge-Kutta methods
-------------------

Explicit methods
~~~~~~~~~~~~~~~~

.. doxygentypedef:: ponio::runge_kutta::euler_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::explicit_euler_sub4_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk6es_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_118_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_21_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_22_midpoint_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_22_ralston_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_32_best_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_33_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_33_233e_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_33_bogackishampine_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_33_heun_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_33_ralston_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_33_van_der_houwen_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_44_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_44_235j_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_44_38_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_44_ralston_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_65_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_65_236a_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_76_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_86_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_nssp_21_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_nssp_32_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_nssp_33_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_nssp_53_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_spp_43_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_ssp_22_heun_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_ssp_32_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_ssp_33_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_ssp_42_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_ssp_53_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk_ssp_54_t
   :project: solver

Embedded methods
~~~~~~~~~~~~~~~~

.. doxygentypedef:: ponio::runge_kutta::rk54_6m_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk54_7m_t
   :project: solver

.. doxygentypedef:: ponio::runge_kutta::rk54_7s_t
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

.. doxygenfunction:: ponio::runge_kutta::leuler
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lexplicit_euler_sub4
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk6es
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_118
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_21
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_22_midpoint
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_22_ralston
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_32_best
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_33
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_33_233e
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_33_bogackishampine
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_33_heun
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_33_ralston
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_33_van_der_houwen
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_44
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_44_235j
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_44_38
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_44_ralston
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_65
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_65_236a
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_76
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_86
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_nssp_21
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_nssp_32
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_nssp_33
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_nssp_53
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_spp_43
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_ssp_22_heun
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_ssp_32
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_ssp_33
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_ssp_42
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_ssp_53
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk_ssp_54
   :project: solver


Embedded methods
~~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::runge_kutta::lrk54_6m
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk54_7m
   :project: solver

.. doxygenfunction:: ponio::runge_kutta::lrk54_7s
   :project: solver

Runge-Kutta Chebyshev methods
-----------------------------

.. doxygenclass:: ponio::runge_kutta::chebyshev::explicit_rkc2
   :project: solver
   :members:
