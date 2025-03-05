Comparison
==========

This section is the comparison with some ODE solvers, we present the same example writing with other ODE solvers to show choices of interface in ponio.

.. toctree::

    ascent/README
    diffeq/README
    gsl/README
    odeint/README
    scipy/README

Lorenz equations
----------------

We would like to solve the Lorenz equations, a classical chaotic system given by:

.. math::

    \begin{aligned}
        \dot{x} &= \sigma ( y - x ) \\
        \dot{y} &= \rho x - y - x z \\
        \dot{z} &= x y - \beta z
    \end{aligned}

with parameter :math:`\sigma = 10`, :math:`\rho = 28` and :math:`\beta = \frac{8}{3}`, and the initial state :math:`(x_0, y_0, z_0) = (1, 1, 1)`. We solve it with a classical Runge-Kutta method of order 4, given by its Butcher tableau:

.. math::

    \begin{array}{c|cccc}
        0           &             &             &             &             \\
        \frac{1}{2} & \frac{1}{2} &             &             &             \\
        \frac{1}{2} & 0           & \frac{1}{2} &             &             \\
        1           & 0           & 0           & 1           &             \\
        \hline
                    & \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6}
    \end{array}

This method is pretty much always in library with various name:

* ``RK4`` in :doc:`Ascent <ascent/README.md>`
* ``RK4`` in :doc:`DifferentialEquations.jl <diffeq/README>`
* ``gsl_odeiv2_step_rk4`` in :doc:`GSL <gsl/README>`
* ``runge_kutta4`` in :doc:`odeint <odeint/README>`
* ``rk_44`` in :doc:`ponio <ponio/README>`

 except in :doc:`SciPy <scipy/README>` where we need to add this.

.. image:: lorenz.png
