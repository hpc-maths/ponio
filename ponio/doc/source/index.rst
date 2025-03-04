.. _topics-index:

=============================
ponio |version| documentation
=============================

.. note::
   ponio means *Presentation Overview of Numerical Integrator for ODE*, and of course it is also a reference to the animated film *Ponyo*.

The library ponio is a collection of time integrators for solving differential equations written in C++. The purpose of this library is to supply efficient and flexible C++ implementation of solvers for various differential equations. Main method classes are :

* explicit Runge-Kutta methods (eRK)
* diagonal implicit Runge-Kutta methods (DIRK)
* Lawson methods (based with an underlying Runge-Kutta method) (LRK)
* exponential Runge-Kutta methods (expRK)
* Runge-Kutta Chebyshev (RKC)
* splitting method (Lie or Strang)

This library aims to be the easiest to use without compromising on performance.

.. toctree::
   :caption: User Documentation
   :maxdepth: 2

   user/user

.. toctree::
   :caption: Developer Documentation
   :maxdepth: 3
   :glob:

   api/api

.. toctree::
   :caption: Gallery
   :glob:

   gallery/gallery
