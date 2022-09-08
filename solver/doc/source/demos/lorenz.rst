Lorenz system
=============

.. math::

  \begin{cases}
    \dot{x} = \sigma(y-x) \\
    \dot{y} = \rho x - y - xz \\
    \dot{z} = xy - \beta z \\
  \end{cases}

with initial condition :math:`(x,y,z)=(1,1,1)` and parameter :math:`\sigma=10`, :math:`\rho=28` and :math:`\beta=\frac{8}{3}`.

Runge-Kutta solver
------------------

First we need to write a problem:

.. code-block:: cpp
  
  using state_t = std::valarray<double>;
  
  state_t lorenz(double t, const state_t& u)
  {
    double sigma=10. , rho=28., beta=8./3.;
    return {
      signa*( u[1] - u[0] ),
      rho*u[0] - u[1] - u[0]*u[2],
      u[0]*u[1] - beta*u[2]
    };
  }

  //...

  auto pb_rk = ode::make_problem(lorenz);

.. note::

  We can also use a lambda to define a problem, or a collection of function as we present in the following part about splitting methods

Next we define the initial condition and call :cpp:func:`ode::solve`:

.. code-block:: cpp

  state_t u0{1.,1.,1.};

  using namespace observer; // to use _fobs litteral

  ode::solve( pb_rk , ode::runge_kutta::rk33() , u0 , {0.,20.} , 0.01 , "lorenz_rk33.dat"_fobs );
    
