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
    ponio/README
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

and the forward Euler method for time discretization:

.. math::

    u^{n+1} = u^n - \Delta t \textrm{D}_a(u^n)

where :math:`\textrm{D}_a` is the approximation of partial derivative in :math:`x` at velocity :math:`a`. The complet scheme, with :math:`a>0`, becomes:

.. math::

    u^{n+1}_i = u^n_i - a \frac{\Delta t}{\Delta t}( u^n_{i+1} - u^n_i )

In case of :math:`a=1` and :math:`\Delta t = \Delta x`, the scheme gives the exact solution. The explicit Euler method (or forward Euler method) is present in:

* :doc:`Ascent <ascent/README>` with the name ``asc::Euler``
* :doc:`odeint <odeint/README>` with the name ``euler``
* :doc:`petsc <petsc/README>` with the name ``TSRK1FE``
* :doc:`ponio <ponio/README>` with the name ``euler``

In :doc:`GSL <gsl/README>` we need to implement our own driver, and in :doc:`SciPy <scipy/README>` to provide a Butcher tableau.

.. image:: transport.png


Arenstorf orbit
---------------

We would like to compare some adaptive time step methods with the Arenstorf orbit problem

.. math::

    \begin{cases}
        \ddot{x} &= x + 2\dot{y} - \frac{1-\mu}{r_1^3}(x+\mu) - \frac{\mu}{r_2^3}(x-1+\mu) \\
        \ddot{y} &= y - 2\dot{x} - \frac{}{1-\mu}{r_1^3}y - \frac{\mu}{r_2^3}y
    \end{cases}

with initial condition :math:`(x,\dot{x},y,\dot{y})=(0.994, 0, 0, -2.001585106)`, :math:`r_1` and :math:`r_2` given by

.. math::

    r_1 = \sqrt{(x+\mu)^2 + y^2},\quad r_2 = \sqrt{(x-1+\mu)^2 + y^2}

and with parameter :math:`\mu = 0.012277471`.

First of all, we need to rewrite this problem into a first order derivative equation in time

.. math::

    \begin{pmatrix}
        y_1 \\
        y_2 \\
        y_3 \\
        y_4
    \end{pmatrix}
    =
    \begin{pmatrix}
        x \\
        y \\
        \dot{x} \\
        \dot{y}
    \end{pmatrix},
    \qquad
    \begin{cases}
        \dot{y}_1 = y_3 \\
        \dot{y}_2 = y_4 \\
        \dot{y}_3 = y_1 + 2y_4 - \frac{1-\mu}{r_1^3}(y_1 + \mu) - \frac{\mu}{r_2^3}(y_1-1+\mu) \\
        \dot{y}_4 = y_2 - 2y_3 - \frac{1-\mu}{r_1^3}y_2 - \frac{\mu}{r_2^3}y_2 \\
    \end{cases}

This is a good example to use an adaptive time step method, because we need a small time step close to the Moon, elsewhere the constraint can be relaxed. For the sake of simplicity we use only native adaptive time step method in each library, it will be:

* ``DOPRI45`` in :doc:`Ascent <ascent/README>`
* ``gsl_odeiv2_step_rk8pd`` in :doc:`GSL <gsl/README>`
* ``runge_kutta_dopri5`` in :doc:`odeint <odeint/README>`
* ``TSRK5DP`` in :doc:`petsc <petsc/README>`
* ``rk54_7m`` or ``rk87_13m`` in :doc:`ponio <ponio/README>`
* ``DOP853`` in :doc:`SciPy <scipy/README>`

where ``DOPRI45``, ``runge_kutta_dopri5`` or ``rk54_7m`` are the same method, in :cite:`dormand:1980` (RK5(4) 7M), see :cpp:type:`ponio::runge_kutta::rk54_7m_t` for more information; and ``gsl_odeiv2_step_rk8pd``, ``DOP853`` and ``rk87_13m`` are the same method in :cite:`prince:1981` (RK8(7) 13M) :cpp:type:`ponio::runge_kutta::rk87_13m_t`.

.. image:: arenstorf.png
