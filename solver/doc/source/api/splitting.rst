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

.. doxygenclass:: ode::splitting::lie
   :project: solver
   :members:

Helper functions and classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: ode::splitting::make_lie_from_tuple
   :project: solver

.. doxygenclass:: ode::splitting::lie_tuple
   :project: solver
   :members:

.. doxygenfunction:: ode::splitting::make_lie_tuple
   :project: solver

Strang splitting method
-----------------------

In Strang splitting method, the solution is computed as:

.. math::

   u^{n+1} = \phi_{\frac{\Delta t}{2}}^{[f_1]}\circ \cdots \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}
              \circ \phi_{\Delta t}^{[f_n]}
              \circ \phi_{\frac{\Delta t}{2}}^{[f_{n-1}]}\circ\cdots\circ \phi_{\frac{\Delta t}{2}}^{[f_1]}
              (t^n,u^n)

.. doxygenclass:: ode::splitting::strang
   :project: solver
   :members:

Helper functions and classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: ode::splitting::make_strang_from_tuple
   :project: solver

.. doxygenclass:: ode::splitting::strang_tuple
   :project: solver
   :members:

.. doxygenfunction:: ode::splitting::make_strang_tuple
   :project: solver
