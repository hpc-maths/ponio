Problems
========

.. digraph:: G
  :caption: UML diagram of problems
  :align: center
  :layout: dot

  compound=true;
  fontname = "Bitstream Vera Sans"
  fontsize = 8

  node [
    fontname = "Bitstream Vera Sans"
    fontsize = 8
    shape = "record"
  ]

  edge [
    fontname = "Bitstream Vera Sans"
    fontsize = 8
  ]

  // classes
  parent_problem [
    label = "{parent_problem}"
  ]
  simple_problem [
    label = "{simple_problem}"
  ]
  lawson_problem [
    label = "{lawson_problem}"
  ]
  imex_problem [
    label = "{imex_problem}"
  ]
  problem [
    label = "{problem}"
  ]

  // inheritence
  edge [ arrowhead = "empty" ]
  simple_problem -> parent_problem
  lawson_problem -> parent_problem
  imex_problem -> parent_problem

  // aggregation
  edge [ arrowhead = odiamond ]
  parent_problem -> problem


Parent problem
--------------

.. doxygenclass:: ode::parent_problem
   :project: solver
   :members:

Simple problem
--------------

.. doxygenclass:: ode::_simple_problem
   :project: solver

.. doxygenclass:: ode::simple_problem
   :project: solver
   :members:

.. doxygenfunction:: ode::make_simple_problem
   :project: solver

Lawson problem
--------------

.. doxygenclass:: ode::_lawson_problem
   :project: solver

.. doxygenclass:: ode::lawson_problem
   :project: solver
   :members:

.. doxygenfunction:: ode::make_lawson_problem
   :project: solver

IMEX problem
------------

.. doxygenclass:: ode::_imex_problem
   :project: solver

.. doxygenclass:: ode::imex_problem
   :project: solver
   :members:

Problem
-------

A problem is a collection of sub-problems defined by previous problems.

.. doxygenclass:: ode::problem
   :project: solver
   :members:
