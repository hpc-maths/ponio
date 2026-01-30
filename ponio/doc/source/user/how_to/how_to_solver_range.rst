How-to: make a :cpp:class:`ponio::solver_range`
----------------------------------------------

The signature of :cpp:func:`ponio::make_solver_range` looks like the same signature of :cpp:func:`ponio::solve` function, without the last parameter (the observer):

* The problem: :code:`f`
* The algorithm used to solve the problem: :cpp:func:`ponio::runge_kutta::rk_33`
* The initial condition: :code:`y0`
* The time interval with a :cpp:class:`ponio::time_span`: :code:`{0., 2.0}`
* The initial time step: :code:`dt`

.. literalinclude:: ../../_static/cpp/how_to_solver_range.cpp
  :language: cpp
  :lines: 24
  :lineno-start: 24
  :linenos:

Next we can get an iterator on :cpp:class:`ponio::solver_range` and loop until the end

.. literalinclude:: ../../_static/cpp/how_to_solver_range.cpp
  :language: cpp
  :lines: 25-31
  :lineno-start: 25
  :linenos:

.. seealso::

  The full example can be found in :download:`how_to_solver_range.cpp <../../_static/cpp/how_to_solver_range.cpp>`.
