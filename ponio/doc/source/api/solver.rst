Solver
======

There are two ways to solve a problem with ponio, the first one is with the function :cpp:func:`ponio::solve`, the second one is with an iterator with :cpp:class:`ponio::time_iterator`.

Function
--------

The function :cpp:func:`ponio::solve` takes care of the time loop, you can't control it. This is sufficient in many cases.

.. doxygenfunction:: ponio::solve
   :project: ponio

Example of solving Curtiss-Hirschfelder problem with RK(4, 4). The ODE is defined by

.. math::

   \begin{cases}
      \dot{y} = k(\cos(t) - y) \\
      y(0) = 2
   \end{cases}

with :math:`k=50`.

.. code-block:: cpp

   auto f = [](double t, double y)
   {
      double k = 50;
      return k * ( std::cos(t) - y );
   };

   double y0 = 2.0;
   double dt = 0.01;

   ponio::solver(f, ponio::runge_kutta::rk_44(), y0, {0., 2.}, dt, "output.txt"_fobs);


----------


Iterator
--------

When you need to control the time loop you have to build a :cpp:struct:`ponio::solver_range` and then iterate on it with a :cpp:struct:`ponio::time_iterator`.

Solver range class
~~~~~~~~~~~~~~~~~~

.. doxygenstruct:: ponio::solver_range
   :project: ponio
   :members:

Helper function for solver_range

.. doxygenfunction:: ponio::make_solver_range
   :project: ponio


Time iterator class
~~~~~~~~~~~~~~~~~~~

.. doxygenstruct:: ponio::time_iterator
   :project: ponio
   :members:


Helper function for time_iterator

.. doxygenfunction:: ponio::make_time_iterator
   :project: ponio

.. doxygenstruct:: ponio::current_solution
   :project: ponio
   :members:


Iteration information class
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need information on current iteration, ponio provides a :cpp:struct:`ponio::iteration_info` class to access to some information on algorithm:

* number of stages,
* number of evaluation of function (also count evaluation in Newton method for implicit methods),
* a boolean to test if iterator is on a step given to :cpp:func:`ponio::make_solver_range` (:cpp:var:`t_span`),
* tolerance, error and a boolean to test if the iteration is succeed for adaptive time step methods.

.. doxygenstruct:: ponio::iteration_info
   :project: ponio
   :members:
