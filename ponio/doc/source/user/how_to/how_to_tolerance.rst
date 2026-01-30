How-to: change tolerances for adaptive time step methods
--------------------------------------------------------

Adaptive time step methods provide an error estimator such as:

.. math::

  \sum_i \left( \frac{|u^{n+1}_i - \tilde{u}^{n+1}_i|}{a_{tol} + \max(|u^n_i|, |u^{n+1}_i||)} \right)

where :math:`u^n`, :math:`u^{n+1}` and :math:`\tilde{u}^{n+1}` are respectively the solution at time :math:`t^n`, time :math:`t^{n+1}=t^n+\Delta t` and an other estimation at time :math:`t^{n+1}=t^n+\Delta t`. We also get two tolerances: the absolute one :math:`a_{tol}` and relative one :math:`r_{tol}`, we can set them with respectively :code:`abs_tol` and :code:`rel_tol` member function of algorithm.

.. literalinclude:: ../../_static/cpp/how_to_tolerance.cpp
  :language: cpp
  :lines: 34-35
  :lineno-start: 34
  :linenos:

You can also change the parameters of Newton method in diagonal implicit method, the tolerance and the maximum number of iteration with respectively :code:`newton_tol` and :code:`newton_max_iter` member functions.

.. literalinclude:: ../../_static/cpp/how_to_tolerance.cpp
  :language: cpp
  :lines: 38-39
  :lineno-start: 38
  :linenos:

.. seealso::

  The full example can be found in :download:`how_to_tolerance.cpp <../../_static/cpp/how_to_tolerance.cpp>`.
