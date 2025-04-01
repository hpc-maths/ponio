Comparison
==========

This section is the comparison with some ODE solvers, we present the same example writing with other ODE solvers to show choices of interface in ponio.

.. toctree::
    :maxdepth: 1

    ascent/README
    diffeq/README
    gsl/README
    odeint/README
    petsc/README
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

* ``RK4`` in :doc:`Ascent <ascent/README>`
* ``RK4`` in :doc:`DifferentialEquations.jl <diffeq/README>`
* ``gsl_odeiv2_step_rk4`` in :doc:`GSL <gsl/README>`
* ``runge_kutta4`` in :doc:`odeint <odeint/README>`
* ``TSRK4`` in :doc:`petsc <petsc/README>`
* ``rk_44`` in :doc:`ponio <ponio/README>`

except in :doc:`SciPy <scipy/README>` where we need to add this.

.. raw:: html
    :file: lorenz.html


Transport equation
------------------

In this example we would like to solve the following PDE:

.. math::

    \partial_t u + a \partial_x u = 0

with :math:`t>0`, on the torus :math:`x\in [0, 1)`, a velocity :math:`a=1` and with the initial condition given by:

.. math::

    u(0, x) = \begin{cases}
        x - 0.25  & \text{if } x\in[0.25, 0.5[ \\
        -x + 0.75 & \text{if } x\in[0.5, 0.75[ \\
        0         & \text{else}
    \end{cases}

We choose a first order up-wind scheme to estimate the :math:`x` derivative:

.. math::

    a \partial_x u(t^n, x_i) \approx \begin{cases}
        \frac{u^n_{i+1} - u^n_{i}}{\Delta x}, & \text{if } a \geq 0 \\
        \frac{u^n_{i} - u^n_{i-1}}{\Delta x}, & \text{else}
    \end{cases}

.. image:: transport.png
