How-to: understand :cpp:class:`ponio::time_iterator`
----------------------------------------------------

After create a :cpp:class:`ponio::solver_range` with :cpp:func:`ponio::make_solver_range` as:

.. literalinclude:: ../../_static/cpp/how_to_solve_range.cpp
  :language: cpp
  :lines: 24
  :lineno-start: 24
  :linenos:

you want to create a :cpp:class:`ponio::time_iterator` to loop over the solver_range

.. literalinclude:: ../../_static/cpp/how_to_solve_range.cpp
  :language: cpp
  :lines: 25-31
  :lineno-start: 25
  :linenos:

The :code:`it_sol` object has the following structure

.. uml:: ../time_iterator.uml
  :caption: Summarize of :code:`time_iterator` class

.. seealso::

  The full example can be found in :download:`how_to_solve_range.cpp <../../_static/cpp/how_to_solve_range.cpp>`.
