Understand observers with Lotka-Volterra equations
==================================================

In this example, we present what can we do with :doc:`observers <../api/observer>` in ponio, with Lotka-Voletrra equations. The equations are the following

.. math::

  \begin{cases}
    \frac{\mathrm{d}x}{\mathrm{d}t} = \alpha x - \beta xy\\
    \frac{\mathrm{d}y}{\mathrm{d}t} = \delta xy - \gamma y\\
  \end{cases}

with :math:`\alpha = \frac{2}{3}`, :math:`\beta = \frac{4}{3}` and :math:`\delta = \gamma = 1`. Initial condition :math:`(x, y)(t=0) = (x_0, x_0)`, and we will take multiple values of :math:`x_0`.

We write the problem like

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.cpp
  :language: cpp
  :lines: 19-31
  :lineno-start: 19
  :linenos:

And launch the simulation between :math:`0` and :math:`15` with the time step :math:`\Delta t = 0.1`

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.cpp
  :language: cpp
  :lines: 36-41
  :lineno-start: 36
  :linenos:


The next sections specify how to define the observer :code:`obs`, to write into a file, in standard output, a stream, or a user-defined observer. The observer is called at each iteration in time and take the current time :math:`t^n`, the current state :math:`u^n` and the current time step :math:`\Delta t^n`.


The file observer
-----------------

In this example we export data into a file with :cpp:class:`observer::file_observer` class. The ponio library provides this output for simple types and containers (thanks concepts). The output is formate as following

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.txt
  :language: text
  :lines: 1-5

First column is the current time, the last one is the current time step, and between we get the multiple values of current state. Other observers provide by ponio have the same format.

The :cpp:class:`observer::file_observer` class can be build with a `std::string <https://en.cppreference.com/w/cpp/string/basic_string>`_ or a `std::filesystem::path <https://en.cppreference.com/w/cpp/filesystem/path>`_. The constructor makes necessary directory if needed.

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.cpp
  :language: cpp
  :lines: 33-34
  :lineno-start: 33
  :linenos:

.. note::

  Data are not flushed into the output file, you have to wait the flush of the buffer.

The full example can be found in :download:`lotka_volterra_fobs.cpp <../_static/cpp/lotka_volterra_fobs.cpp>`.



The :code:`cout` observer
-------------------------

In this example we export data into the standard output with :cpp:class:`observer::cout_observer` class.

.. literalinclude:: ../_static/cpp/lotka_volterra_cobs.cpp
  :language: cpp
  :lines: 31
  :lineno-start: 31
  :linenos:

.. note::

  Data are not flushed into the standard output, you have to wait the flush of the buffer.

The full example can be found in :download:`lotka_volterra_cobs.cpp <../_static/cpp/lotka_volterra_cobs.cpp>`.



The stream observer
-------------------

In this example we export data into the standard output with :cpp:class:`observer::stream_observer` class. This observer is the more generally one and can be plug with any other type of stream, for example a `std::stringstream <https://en.cppreference.com/w/cpp/io/basic_stringstream>`_. The :cpp:class:`observer::file_observer` is a specialization with a `std::ofstream <https://en.cppreference.com/w/cpp/io/basic_ofstream>`_ and :cpp:class:`observer::cout_observer` with `std::cout <https://en.cppreference.com/w/cpp/io/cout>`_.


.. literalinclude:: ../_static/cpp/lotka_volterra_sobs.cpp
  :language: cpp
  :lines: 33-34
  :lineno-start: 33
  :linenos:

Next we get all informations into the observer or our buffer.

.. literalinclude:: ../_static/cpp/lotka_volterra_sobs.cpp
  :language: cpp
  :lines: 43
  :lineno-start: 43
  :linenos:


The full example can be found in :download:`lotka_volterra_sobs.cpp <../_static/cpp/lotka_volterra_sobs.cpp>`.



The user-defined observer
-------------------------

An observer is a invocable object (function, functor, lambda) with the following signature

.. code-block:: cpp

  void operator() ( value_t current_time, state_t & current_state, value_t current_time_step );

In our case, :code:`value_t` is :code:`double` and :code:`state_t` is :code:`std::valarray<double>`. We propose to export also the invariant of Lotka-Voletrra equations given by

.. math::

  V = \delta x - \ln(x) + \beta y - \alpha \ln(y)

We implement the observer into a class

.. literalinclude:: ../_static/cpp/lotka_volterra_uobs.cpp
  :language: cpp
  :lines: 13-50
  :lineno-start: 13
  :linenos:

Now we just have to build an instance of this class

.. literalinclude:: ../_static/cpp/lotka_volterra_uobs.cpp
  :language: cpp
  :lines: 71
  :lineno-start: 71
  :linenos:

The full example can be found in :download:`lotka_volterra_uobs.cpp <../_static/cpp/lotka_volterra_uobs.cpp>`.
