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
