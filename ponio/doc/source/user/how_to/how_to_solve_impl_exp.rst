How-to: solve an ODE with an implicit method
--------------------------------------------

.. note::

  You can also see `Very first steps with Lorenz system <../lorenz>`_ tutorial.

We would like to solve the following equation

.. math::

  \dot{y} = -y,\qquad y(0) = 1

We first define the function :math:`f:(t,y)\mapsto f(t, y) = -y`, in ponio we define this function as :code:`f(double t, double y, double& dy)`. We also need to define the Jacobian function of this function :math:`J_f:(t,y)\mapsto \partial_y f(t,y) = -1`, in ponio we define this function as :code:`df(double t, double y)->double`. We store this two functions in a :cpp:class:`ponio::implicit_problem` class:

.. literalinclude:: ../../_static/cpp/how_to_solve_impl_exp.cpp
  :language: cpp
  :lines: 20-28
  :lineno-start: 20
  :linenos:

We also define initial condition and time step:

.. literalinclude:: ../../_static/cpp/how_to_solve_impl_exp.cpp
  :language: cpp
  :lines: 30-31
  :lineno-start: 30
  :linenos:

Next we solve this equation with explicit Euler method defines in ponio as :cpp:type:`ponio::runge_kutta::euler`.

.. literalinclude:: ../../_static/cpp/how_to_solve_impl_exp.cpp
  :language: cpp
  :lines: 33
  :lineno-start: 33
  :linenos:

In this line we call :cpp:func:`ponio::solve` function, with a :cpp:class:`observer::file_observer` to store the output in :code:`how_to_solve_impl_exp.txt` text file.

.. seealso::

  The full example can be found in :download:`how_to_solve_impl_exp.cpp <../../_static/cpp/how_to_solve_impl_exp.cpp>`.
