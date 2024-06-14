Splitting methods
=================

Some methods to solve a problem of the form:

.. math::

  \dot{u} = \sum_i f_i(t,u)

with the initial condition :math:`u(t=0)=u_0`. We note with :math:`\phi_{\tau}^{[f_i]}(t^n,\tilde{u}^n)` the solution at time :math:`t^n+\tau` of the subproblem :math:`\dot{u}=f_i(t,u)` with the initial condition :math:`u(t^n)=\tilde{u}^n`.

Lie splitting method
--------------------

In Lie splitting method, the solution is computed as:

.. math::

   u^{n+1} = \phi_{\Delta t}^{[f_1]}\circ \cdots \circ \phi_{\Delta t}^{[f_n]} (t^n,u^n)

.. doxygenclass:: ponio::splitting::lie::lie
   :project: ponio
   :members:

Helper function
~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::splitting::lie::make_lie_tuple
   :project: ponio

Strang splitting method
-----------------------

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}
              \circ \phi_{\Delta t}^{[f_n]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

.. doxygenclass:: ponio::splitting::strang::strang
   :project: ponio
   :members:

Helper function
~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::splitting::strang::make_strang_tuple
   :project: ponio


Adaptive time step Strang splitting method
------------------------------------------

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}
              \circ \phi_{\Delta t}^{[f_n]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

The approximation of order 1 is computed with a shifted Strang splitting method:

.. math::

   \tilde{u}^{n+1} = \phi_{(\frac{1}{2}-\delta)\Delta t}^{[f_1]}\circ\phi_{\frac{\Delta t}{2}}^{[f_2]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}
              \circ \phi_{\Delta t}^{[f_n]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}\circ\cdots\circ\phi_{\frac{\Delta t}{2}}^{[f_2]}\circ \phi_{(\frac{1}{2}+\delta)\Delta t}^{[f_1]}
              (t^n,u^n)

with :math:`\delta\in[-1/2, 0)\cup(0,1/2]`

.. note::

   Only the first (and last) step is shifted by :math:`\delta\Delta t`.

.. doxygenclass:: ponio::splitting::strang::adaptive_strang
   :project: ponio
   :members:

Helper function
~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::splitting::strang::make_adaptive_strang_tuple
   :project: ponio

Generic splitting tuple
-----------------------

The class :cpp:class:`ponio::splitting::detail::_splitting_tuple` is useful to build each splitting method from variadic number of pair of algorithm and time step. This class can has optional arguments, stored in a tuple, to call constructor of :code:`_splitting_method_t` (which is :cpp:class:`ponio::splitting::lie::lie`, :cpp:class:`ponio::splitting::strang::strang` or :cpp:class:`ponio::splitting::strang::adaptive_strang`).

.. doxygenclass:: ponio::splitting::detail::_splitting_tuple
   :project: ponio
   :members:

Helper function
~~~~~~~~~~~~~~~

.. doxygenfunction:: ponio::splitting::detail::make_splitting_from_tuple( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts )
   :project: ponio

.. doxygenfunction:: ponio::splitting::detail::make_splitting_from_tuple( std::tuple<Methods_t...> const& meths, std::array<value_t, sizeof...( Methods_t )> const& dts, optional_tuple_t optional_args )
   :project: ponio
