Method class
============

.. note::
  **Naming convention:** about type name of template parameter.

  * ``Algorithm_t`` matches to an algorithm *i.e.* a class that contains stages of a Runge-Kutta method. It needs to provide an operator ``()`` with the following definition:

    .. code-block:: cpp

      state_t operator()( Stage<s>, Problem_t& f, value_t tn, state_t const& un, ArrayKi_t const& ki, value_t dt)
  * ``Method_t`` matches to the ``ode::method`` class, and provides the following operator ``()``:

    .. code-block:: cpp

      state_t operator()(Problem_t& f, value_t tn, state_t const& un, value_t dt)

.. doxygenclass:: ode::method
   :project: solver
   :members:

There is 0 overloading of ``ode::make_method`` function, one for generic algorithm, and one for each implemanted splitting method (Lie splitting method and Strang splitting method).

.. doxygenfunction:: ode::make_method(Algorithm_t const& algo, state_t const& shadow_of_u0)
   :project: solver

.. doxygenfunction:: ode::make_method(splitting::lie_tuple<Algorithms_t...> const &algos, state_t const &shadow_of_u0)
   :project: solver

.. doxygenfunction:: ode::make_method(splitting::strang_tuple<Algorithms_t...> const& algos, state_t const& shadow_of_u0)
   :project: solver

An helper function for splitting method factory.

.. doxygenfunction:: ode::make_tuple_methods
   :project: solver

