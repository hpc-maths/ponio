Very first steps with Lorenz system
===================================

In this example, we present how to solve Lorenz system, step by step, with ponio library. The Lorenz system is the following

.. math::

  \begin{cases}
    \dot{x} = \sigma(y-x) \\
    \dot{y} = \rho x - y - xz \\
    \dot{z} = xy - \beta z \\
  \end{cases}

with initial condition :math:`(x,y,z)=(1,1,1)` and parameter :math:`\sigma=10`, :math:`\rho=28` and :math:`\beta=\frac{8}{3}`. We present in this section how to solve it with an explicit Runge-Kutta method, a Lawson method, a diagonal implicit Runge-Kutta method and a splitting method.


Explicit Runge-Kutta solver
---------------------------

In this section we solve this problem with a classical Runge-Kutta method. First of all we need to rewrite the problem as :math:`\dot{u} = f(t, u)`, here we choose :math:`u = (x, y, z)`

.. math::

  u = \begin{pmatrix}
    x \\
    y \\
    z \\
  \end{pmatrix}, \quad
  \dot{u} = f(t, u) = \begin{pmatrix}
    \sigma(u_1 - u_0) \\
    \rho u_0 - u_1 - u_0 u_2 \\
    u_0 u_1 - \beta u_2
  \end{pmatrix}

next we choose a data structure to store state :math:`u`, for the sake of simplicity we chose the STL container :code:`std::valarray<double>`:

.. literalinclude:: ../_static/cpp/lorenz_rk.cpp
  :language: cpp
  :lines: 16
  :lineno-start: 16
  :linenos:

Now we can write :math:`f:t, u\mapsto f(t,u)` as a lambda function:

.. literalinclude:: ../_static/cpp/lorenz_rk.cpp
  :language: cpp
  :lines: 18-27
  :lineno-start: 18
  :linenos:

we can also use a :cpp:class:`ponio::problem` object, a user-defined functor or a function. Next we define the initial condition to :math:`u_0 = (1, 1, 1)`, the initial and final time in a time span (here between :math:`0` and :math:`20`) the time step :math:`\Delta t = 0.01` and solve with all this parameters with a classical RK(4,4): :cpp:type:`ponio::runge_kutta::rk_44_t`.

.. note::

  List of all explicit Runge-Kutta methods can be found in the :doc:`algorithms section <../api/algorithm>`.

.. literalinclude:: ../_static/cpp/lorenz_rk.cpp
  :language: cpp
  :lines: 29-34
  :lineno-start: 29
  :linenos:
  :emphasize-lines: 6

In the last line we call :cpp:func:`ponio::solve` function, with a :cpp:class:`observer::file_observer` to store the output in :code:`lorenz_rk44.txt` text file.

.. note::

  To use the :code:`_fobs` literal you need to using the correct namespace

  .. literalinclude:: ../_static/cpp/lorenz_rk.cpp
    :language: cpp
    :lines: 15

The full example can be found in :download:`lorenz_rk.cpp <../_static/cpp/lorenz_rk.cpp>`.

You can compile this example

.. code-block::

  $CXX -std=c++20 -I PONIO_INCLUDE lorenz_rk.cpp -o lorenz_rk

or with a cmake like in `ponio-gallery <https://github.com/hpc-maths/ponio-gallery>`_ repository.


----


Lawson solver
-------------

In this section we propose to solve the Lorenz equations with a Lawson method. First of all we need to rewrite the system with a linear and non-linear part, as :math:`\dot{u} = Lu + N(t, u)`, with :math:`u=(x, y, z)`

.. math::

  u = \begin{pmatrix}
    x \\
    y \\
    z \\
  \end{pmatrix}, \quad
  \dot{u} = \underbrace{\begin{pmatrix}
    -\sigma & \sigma & 0 \\
    \rho    & -1     & 0 \\
    0       & 0      & -\beta
  \end{pmatrix}}_{L}
  u
  +
  \underbrace{\begin{pmatrix}
    0 \\
    -u_0 u_2 \\
    u_0 u_1 \\
  \end{pmatrix}}_{N(t,u)}

A Lawson method needs to compute an exponential of the linear part, but the ponio library doesn't provide linear algebra, and even less an exponential function for matrices. For that you can use `Eigen 3 <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_, so the definition of type of state :math:`u^n` changes into

.. literalinclude:: ../_static/cpp/lorenz_lrk.cpp
  :language: cpp
  :lines: 18-21
  :lineno-start: 18
  :linenos:

Now we can write the linear part :math:`L` and non-linear part :math:`N:t, u\mapsto N(t,u)`, and store them into a :cpp:class:`ponio::lawson_problem`

.. literalinclude:: ../_static/cpp/lorenz_lrk.cpp
  :language: cpp
  :lines: 23-34
  :lineno-start: 23
  :linenos:
  :emphasize-lines: 12

The :cpp:class:`ponio::lawson_problem` builds in the emphasize line provides an interface which can be use by a classical explicit Runge-Kutta method or a Lawson method.

We have to define the exponential function which takes a matrix and return a matrix

.. literalinclude:: ../_static/cpp/lorenz_lrk.cpp
  :language: cpp
  :lines: 36-40
  :lineno-start: 36
  :linenos:

We can also define your own exponential function if you know the analytic form of the result of :math:`\exp(\tau L)` with :math:`\tau = c_j\Delta t`.

The call of function :cpp:func:`ponio::solve` doesn't change

.. literalinclude:: ../_static/cpp/lorenz_lrk.cpp
  :language: cpp
  :lines: 42-47
  :lineno-start: 42
  :linenos:

We use the Lawson method given by the underlying RK(4,4) method: :cpp:var:`ponio::runge_kutta::lrk_44_t`.

.. note::

  List of all Lawson methods can be found in the :doc:`algorithms section <../api/algorithm>`.

The full example can be found in :download:`lorenz_lrk.cpp <../_static/cpp/lorenz_lrk.cpp>`.


----


Diagonal implicit Runge-Kutta solver
------------------------------------

In this section we solve this problem with an implicit Runge-Kutta method. We keep the problem with the form :math:`\dot{u} = f(t,u)` but we need to provide the Jacobian function for the Newton method at each stage.

Because we will use linear algebra, we use `Eigen 3 <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_ to define the state

.. literalinclude:: ../_static/cpp/lorenz_dirk.cpp
  :language: cpp
  :lines: 19-22
  :lineno-start: 19
  :linenos:

Now we define the function and its Jacobian defined by

.. math::

  \partial_u f = \begin{pmatrix}
    -\sigma    & \sigma & 0    \\
    \rho - u_2 & - 1    & -u_0 \\
    u_1        & u_0    & -\beta
  \end{pmatrix}

and create a :cpp:class:`ponio::implicit_problem`

.. literalinclude:: ../_static/cpp/lorenz_dirk.cpp
  :language: cpp
  :lines: 24-42
  :lineno-start: 24
  :linenos:

Finally we can solve the problem with the same lines of previous solver

.. literalinclude:: ../_static/cpp/lorenz_dirk.cpp
  :language: cpp
  :lines: 44-49
  :lineno-start: 44
  :linenos:

We use the DIRK(3, 4) method: :cpp:func:`ponio::runge_kutta::dirk34_t`.

.. note::

  List of all DIRK methods can be found in the :doc:`algorithms section <../api/algorithm>`.

The full example can be found in :download:`lorenz_dirk.cpp <../_static/cpp/lorenz_dirk.cpp>`.


----


Splitting solver
----------------

For the splitting method, we split the problem into multiple parts. We chose here 3 parts, but ponio can take any number of part for the splitting. We choose arbitrarily the following splitting

.. math::

  \varphi^{[0]} = \begin{cases}
    \dot{x} &= \sigma y \\
    \dot{y} &= \rho x \\
    \dot{z} &= xy
  \end{cases}\ ,
  \qquad
  \varphi^{[1]} = \begin{cases}
    \dot{x} &= -\sigma x \\
    \dot{y} &= -y \\
    \dot{z} &= -\beta z
  \end{cases}\ ,
  \qquad
  \varphi^{[2]} = \begin{cases}
    \dot{x} &= 0 \\
    \dot{y} &= -xz \\
    \dot{z} &= 0
  \end{cases}

The splitting is defined as

.. math::

  u^{n+1} = \varphi_{\Delta t}^{[0]} \circ \varphi_{\Delta t}^{[1]} \circ \varphi_{\Delta t}^{[2]}(u^n)

for the Lie splitting method, and

.. math::

  u^{n+1} = \varphi_{\Delta t/2}^{[0]} \circ \varphi_{\Delta t/2}^{[1]} \circ \varphi_{\Delta t}^{[2]} \circ \varphi_{\Delta t/2}^{[1]} \circ \varphi_{\Delta t/2}^{[0]}(u^n)

for the Strang splitting method, where :math:`\varphi_{\tau}^{[j]}` is a numerical approximation of the step :math:`j` between :math:`0` and :math:`\tau` with a specific time step :math:`\delta t < \tau`. For more information see :doc:`splitting section <../api/splitting>`.

For the sake of simplicity, we solve each step with a explicit Runge-Kutta method, but you can use any other solver (even a splitting method!).

Like other solver, first of all we define our state

.. literalinclude:: ../_static/cpp/lorenz_split.cpp
  :language: cpp
  :lines: 18
  :lineno-start: 18
  :linenos:

Next we define the three function to represent :math:`\varphi^{[0]}`, :math:`\varphi^{[1]}`, :math:`\varphi^{[2]}`

.. literalinclude:: ../_static/cpp/lorenz_split.cpp
  :language: cpp
  :lines: 20-46
  :lineno-start: 20
  :linenos:

Now we create a tuple of algorithms with :cpp:func:`ponio::splitting::lie::make_lie_tuple` or :cpp:func:`ponio::splitting::strang::make_strang_tuple` for respectively a Lie splitting method or Strang splitting method.

.. literalinclude:: ../_static/cpp/lorenz_split.cpp
  :language: cpp
  :lines: 48-50
  :lineno-start: 48
  :linenos:

We have to select an algorithm for each step of splitting method, and specify a time step (smaller than splitting time step), in this case we choose 3 different Runge-Kutta methods of order 4 with 4 stages:

* A first RK(4,4) 3/8 rule: :cpp:type:`ponio::runge_kutta::rk_44_38_t`, with a time step :math:`\delta t = 0.01` (equals to the splitting time step)
* A second RK(4,4), classical one: :cpp:type:`ponio::runge_kutta::rk_44_t`, with a time step :math:`\delta t = 0.005`
* A third RK(4,4), from Ralston: :cpp:type:`ponio::runge_kutta::rk_44_ralston_t`, with a time step :math:`\delta t = 0.0005`

For each algorithm we specify a time step to solve each :math:`\varphi^{[i]}`.

For a Lie splitting this code become

.. code-block:: cpp

  auto lie = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_44_38(), 0.01 ),
        std::make_pair( ponio::runge_kutta::rk_44(), 0.005 ),
        std::make_pair( ponio::runge_kutta::rk_44_ralston(), 0.0005 ) );

Finally we call :cpp:func:`ponio::solve` function, as previous solver

.. literalinclude:: ../_static/cpp/lorenz_split.cpp
  :language: cpp
  :lines: 52-57
  :lineno-start: 52
  :linenos:


The full example can be found in :download:`lorenz_split.cpp <../_static/cpp/lorenz_split.cpp>`.
