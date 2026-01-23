Before starting
===============

Installation
------------

**From conda:**

.. code-block:: sh

  conda install conda-forge::ponio


**From sources:**

.. code-block:: sh

  git clone https://github.com/hpc-maths/ponio.git
  cd ponio

and install dependencies from conda

.. code-block:: sh

  conda env create -f environment/conda-environment.yml

next install ponio

.. code-block:: sh

  cmake . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build --target install

.. note::

  You can install and launch all examples and theirs dependencies with following lines

  .. code-block:: sh

    conda env create -f environment/conda-environment-all.yml

  and compile all examples with

  .. code-block:: sh

    cmake . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_ALL_EXAMPLES=ON
    cd build
    make visu # to launch all examples and script to plot output

.. tip::

  For every example in ponio their is a source file in C++ and a Python script to compile, launch and plot results. You can compile an example with

  .. code-block:: sh

    make lotka_volterra

  You can also launch the Python script associated with the example to see results with

  .. code-block:: sh

    make visu_lotka_volterra
    # or
    make lotka_volterra_visu

You can also install dependencies with `pixi <https://pixi.prefix.dev/latest/>`_ with

.. code-block:: sh

  pixi install


List of variables available in ``CMakeLists.txt``
-------------------------------------------------

========================== =============  =====================================
Variable                   Default value  Description
========================== =============  =====================================
``BUILD_TESTS``            ``OFF``        Set to ``ON`` to build the tests with `doctest <https://github.com/doctest/doctest>`_
``BUILD_EXAMPLES``         ``OFF``        Set to ``ON`` to build examples (without other dependencies)
``BUILD_EIGEN_EXAMPLES``   ``OFF``        Set to ``ON`` to build examples with `Eigen <https://libeigen.gitlab.io/>`_
``BUILD_CLI11_EXAMPLES``   ``OFF``        Set to ``ON`` to build examples with `CLI11 <https://github.com/CLIUtils/CLI11>`_
``BUILD_SAMURAI_EXAMPLES`` ``OFF``        Set to ``ON`` to build examples with `samurai <https://github.com/hpc-maths/samurai>`_
``BUILD_ALL_EXAMPLES``     ``OFF``        Set to ``ON`` to build all examples with extra dependencies
========================== =============  =====================================

``BUILD_TESTS``
  Set to ``ON`` to build the tests with `doctest <https://github.com/doctest/doctest>`_ (default ``OFF``).

``BUILD_EXAMPLES``
  Set to ``ON`` to build examples (without other dependencies) (default ``OFF``).

``BUILD_EIGEN_EXAMPLES``
  Set to ``ON`` to build examples with `Eigen <https://libeigen.gitlab.io/>`_ (default ``OFF``).

``BUILD_CLI11_EXAMPLES``
  Set to ``ON`` to build examples with `CLI11 <https://github.com/CLIUtils/CLI11>`_ (default ``OFF``).

``BUILD_SAMURAI_EXAMPLES``
  Set to ``ON`` to build examples with `samurai <https://github.com/hpc-maths/samurai>`_ (default ``OFF``).

``BUILD_ALL_EXAMPLES``
  Set to ``ON`` to build all examples with extra dependencies (default ``OFF``).
