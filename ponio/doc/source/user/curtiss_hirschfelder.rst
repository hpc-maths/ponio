Manage your time loop with Curtiss-Hirschfelder equation
========================================================

In this example we solve the Curtiss-Hirschfelder equation

.. math::

  \dot{y} = k(\cos(t) - y)

with stiff parameter :math:`k=50` and initial condition :math:`y(0) = y_0 = 2`. As the name of the :math:`k` parameter suggests, this is a stiff equation, so the time step could be constraint by the value of this parameter if we solve this equation with an explicit Runge-Kutta method.

In this case, we need a small time step at the beginning of the simulation, but next, close to the equilibrium we can get larger time step. We can't do that with :cpp:func:`ponio::solve` function because we need to change the current step during the simulation. We will use a :cpp:class:`ponio::solver_range` class and access to it with a :cpp:class:`ponio::time_iterator` to manage the time loop.

First of all, we write the problem, in this case the state will be stored into a :code:`double`, so we don't specify a :code:`state_t` for the sake of simplicity.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_solve.cpp
  :language: cpp
  :lines: 15-20
  :lineno-start: 15
  :linenos:

Next we prepare our simulation, time span, time step, and initial condition

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_solve.cpp
  :language: cpp
  :lines: 22-25
  :lineno-start: 22
  :linenos:


:code:`ponio::solve` function
-----------------------------

The standard function to solve a ODE in ponio is with the :cpp:func:`ponio::solve` function

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_solve.cpp
  :language: cpp
  :lines: 27-29
  :lineno-start: 27
  :linenos:

The full example can be found in :download:`curtiss_hirschfelder_solve.cpp <../_static/cpp/curtiss_hirschfelder_solve.cpp>`.

.. figure:: ../_static/cpp/curtiss_hirschfelder_solve.png
    :width: 500 px
    :alt: Solution and time step history

    Solution and time step history


A while loop
------------

With the :cpp:func:`ponio::solve` function we can't manage solution, or time step. To do this we need a :cpp:class:`ponio::solver_range`

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_while.cpp
  :language: cpp
  :lines: 29
  :lineno-start: 29
  :linenos:

This object is range in C++20 meaning, so we can iterate on it with a :cpp:class:`ponio::time_iterator` object

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_while.cpp
  :language: cpp
  :lines: 30
  :lineno-start: 30
  :linenos:

This kind of iterator can be incremented (:code:`++it_sol`), and we can access to its stored data with :code:`*` operator (:code:`*it_sol`) and get a tuple with :math:`(t^n, y^n, \Delta t^n)`. We can also access with :code:`->` operator with :code:`it_sol->time` for :math:`t^n`, :code:`it_sol->state` for :math:`y^n` and :code:`it_sol->time_step` for :math:`\Delta t^n`.

To manage the time loop we write a while loop

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_while.cpp
  :language: cpp
  :lines: 32-47
  :lineno-start: 32
  :linenos:

We can access in read and write to data in iterator.

The full example can be found in :download:`curtiss_hirschfelder_while.cpp <../_static/cpp/curtiss_hirschfelder_while.cpp>`.

.. figure:: ../_static/cpp/curtiss_hirschfelder_while.png
    :width: 500 px
    :alt: Solution and time step history

    Solution and time step history

A for loop
----------

Like for a while loop, we need to build a :cpp:class:`ponio::solver_range` and iterate on it by using a range-based for loop.

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_for.cpp
  :language: cpp
  :lines: 29
  :lineno-start: 29
  :linenos:

If you want to modify the tuple :math:`(t^n, y^n, \Delta t^n)` we need a reference :

.. code-block:: cpp

  for ( auto& ui : sol_range )
  // ...

Else we can just use a value

.. code-block:: cpp

  for ( auto ui : sol_range )
  // ...

The complet for loop is the following

.. literalinclude:: ../_static/cpp/curtiss_hirschfelder_for.cpp
  :language: cpp
  :lines: 31-43
  :lineno-start: 31
  :linenos:

The full example can be found in :download:`curtiss_hirschfelder_for.cpp <../_static/cpp/curtiss_hirschfelder_for.cpp>`.

.. figure:: ../_static/cpp/curtiss_hirschfelder_for.png
    :width: 500 px
    :alt: Solution and time step history

    Solution and time step history
