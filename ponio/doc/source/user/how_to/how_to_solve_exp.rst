How-to: solve a simple ODE
--------------------------

.. note::

  You can also see `Very first steps with Lorenz system <../lorenz>`_ tutorial.

We would like to solve the following equation

.. math::

  \dot{y} = -y,\qquad y(0) = 1

We first define the function :math:`f:(t,y)\mapsto f(t, y) = -y`, in ponio we define this function as :code:`f(double t, double y, double& dy)`, where the current state :code:`y` can be capture by value, reference or constant reference.

.. literalinclude:: ../../_static/cpp/how_to_solve_exp.cpp
  :language: cpp
  :lines: 19-22
  :lineno-start: 19
  :linenos:

We also define initial condition and time step:

.. literalinclude:: ../../_static/cpp/how_to_solve_exp.cpp
  :language: cpp
  :lines: 24-25
  :lineno-start: 24
  :linenos:


Next we solve this equation with explicit Euler method defines in ponio as :cpp:type:`ponio::runge_kutta::euler`.

.. literalinclude:: ../../_static/cpp/how_to_solve_exp.cpp
  :language: cpp
  :lines: 27
  :lineno-start: 27
  :linenos:

In this line we call :cpp:func:`ponio::solve` function, with a :cpp:class:`observer::file_observer` to store the output in :code:`how_to_solve_exp.txt` text file.

.. note::

  To use the :code:`_fobs` literal you need to using the correct namespace

  .. literalinclude:: ../../_static/cpp/how_to_solve_exp.cpp
    :language: cpp
    :lines: 17

.. seealso::

  The full example can be found in :download:`how_to_solve_exp.cpp <../../_static/cpp/how_to_solve_exp.cpp>`.
