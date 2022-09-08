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

This Lawson methods have a underlying Runge-Kutta method, they have the same name just the namespace change.

Base class
~~~~~~~~~~

.. doxygenstruct:: ode::lawson::base_rk
   :project: solver
   :members:

