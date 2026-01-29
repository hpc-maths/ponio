How-to guides
=============


.. toctree::
   :caption: How-to guides

   how_to


How-to: solve a simple ODE
--------------------------

.. note::

  You can also see `Very first steps with Lorenz system <lorenz>`_ tutorial.

We would like to solve the following equation

.. math::

  \dot{y} = -y,\qquad y(0) = 1

We first define the function :math:`f:(t,y)\mapsto f(t, y) = -y`, in ponio we define this function as :cpp:`f(double t, double y, double& dy)`, where the current state :cpp:`y` can be capture by value, reference or constant reference.

.. literalinclude:: ../_static/cpp/how_to_solve_exp.cpp
  :language: cpp
  :lines: 19-22
  :lineno-start: 19
  :linenos:

We also define initial condition and time step:

.. literalinclude:: ../_static/cpp/how_to_solve_exp.cpp
  :language: cpp
  :lines: 24-25
  :lineno-start: 24
  :linenos:


Next we solve this equation with explicit Euler method defines in ponio as :cpp:type:`ponio::runge_kutta::euler`.

.. literalinclude:: ../_static/cpp/how_to_solve_exp.cpp
  :language: cpp
  :lines: 27
  :lineno-start: 27
  :linenos:

In this line we call :cpp:func:`ponio::solve` function, with a :cpp:class:`observer::file_observer` to store the output in :code:`how_to_solve_exp.txt` text file.

.. note::

  To use the :code:`_fobs` literal you need to using the correct namespace

  .. literalinclude:: ../_static/cpp/how_to_solve_exp.cpp
    :language: cpp
    :lines: 17

.. seealso::

  The full example can be found in :download:`how_to_solve_exp.cpp <../_static/cpp/how_to_solve_exp.cpp>`.

---

How-to: solve an ODE with an implicit method
--------------------------------------------

.. note::

  You can also see `Very first steps with Lorenz system <lorenz>`_ tutorial.

We would like to solve the following equation

.. math::

  \dot{y} = -y,\qquad y(0) = 1

We first define the function :math:`f:(t,y)\mapsto f(t, y) = -y`, in ponio we define this function as :code:`f(double t, double y, double& dy)`. We also need to define the Jacobian function of this function :math:`J_f:(t,y)\mapsto \partial_y f(t,y) = -1`, in ponio we define this function as :code:`df(double t, double y)->double`. We store this two functions in a :cpp:class:`ponio::implicit_problem` class:

.. literalinclude:: ../_static/cpp/how_to_solve_impl_exp.cpp
  :language: cpp
  :lines: 20-28
  :lineno-start: 20
  :linenos:

We also define initial condition and time step:

.. literalinclude:: ../_static/cpp/how_to_solve_impl_exp.cpp
  :language: cpp
  :lines: 30-31
  :lineno-start: 30
  :linenos:

Next we solve this equation with explicit Euler method defines in ponio as :cpp:type:`ponio::runge_kutta::euler`.

.. literalinclude:: ../_static/cpp/how_to_solve_impl_exp.cpp
  :language: cpp
  :lines: 33
  :lineno-start: 33
  :linenos:

In this line we call :cpp:func:`ponio::solve` function, with a :cpp:class:`observer::file_observer` to store the output in :code:`how_to_solve_impl_exp.txt` text file.

.. seealso::

  The full example can be found in :download:`how_to_solve_impl_exp.cpp <../_static/cpp/how_to_solve_impl_exp.cpp>`.


How-to: make a :cpp:class:`ponio::solver_range`
----------------------------------------------

This looks like the same signature of :cpp:func:`ponio::solve` function, without the last parameter (the observer):

* The problem
* The algorithm used to solve the problem
* The initial condition
* The time interval
* The initial time step

.. literalinclude:: ../_static/cpp/how_to_solve_range.cpp
  :language: cpp
  :lines: 24
  :lineno-start: 24
  :linenos:

Next we can get an iterator on :cpp:class:`ponio::solver_range` and loop until the end

.. literalinclude:: ../_static/cpp/how_to_solve_range.cpp
  :language: cpp
  :lines: 25-31
  :lineno-start: 25
  :linenos:

.. seealso::

  The full example can be found in :download:`how_to_solve_range.cpp <../_static/cpp/how_to_solve_range.cpp>`.


How-to: change tolerances for adaptive time step methods
--------------------------------------------------------

Adaptive time step methods provide an error estimator such as:

.. math::

  \sum_i \left( \frac{|u^{n+1}_i - \tilde{u}^{n+1}_i|}{a_{tol} + \max(|u^n_i|, |u^{n+1}_i||)} \right)

where :math:`u^n`, :math:`u^{n+1}` and :math:`\tilde{u}^{n+1}` are respectively the solution at time :math:`t^n`, time :math:`t^{n+1}=t^n+\Delta t` and an other estimation at time :math:`t^{n+1}=t^n+\Delta t`. We also get two tolerances: the absolute one :math:`a_{tol}` and relative one :math:`r_{tol}`, we can set them with respectively :code:`abs_tol` and :code:`rel_tol` member function of algorithm.

.. literalinclude:: ../_static/cpp/how_to_tolerance.cpp
  :language: cpp
  :lines: 34-35
  :lineno-start: 34
  :linenos:

You can also change the parameters of Newton method in diagonal implicit method, the tolerance and the maximum number of iteration with respectively :code:`newton_tol` and :code:`newton_max_iter` member functions.

.. literalinclude:: ../_static/cpp/how_to_tolerance.cpp
  :language: cpp
  :lines: 38-39
  :lineno-start: 38
  :linenos:

.. seealso::

  The full example can be found in :download:`how_to_tolerance.cpp <../_static/cpp/how_to_tolerance.cpp>`.
