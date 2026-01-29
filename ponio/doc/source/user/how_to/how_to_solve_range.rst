How-to: make a :cpp:class:`ponio::solver_range`
----------------------------------------------

This looks like the same signature of :cpp:func:`ponio::solve` function, without the last parameter (the observer):

* The problem
* The algorithm used to solve the problem
* The initial condition
* The time interval
* The initial time step

.. literalinclude:: ../../_static/cpp/how_to_solve_range.cpp
  :language: cpp
  :lines: 24
  :lineno-start: 24
  :linenos:

Next we can get an iterator on :cpp:class:`ponio::solver_range` and loop until the end

.. literalinclude:: ../../_static/cpp/how_to_solve_range.cpp
  :language: cpp
  :lines: 25-31
  :lineno-start: 25
  :linenos:

.. seealso::

  The full example can be found in :download:`how_to_solve_range.cpp <../../_static/cpp/how_to_solve_range.cpp>`.
