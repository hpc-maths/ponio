# README

This notebooks are writing to explain how to start with Ponio. To compile and launch simulation with Ponio we need some dependancies:

* a C++ 20 compiler
* Jupyter with a complet scientific Python environment:
	+ NumPy
	+ Matplotlib
	+ SciPy
	+ SymPy

Some examples needs [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).

```sh
conda create env ponio-notebook
conda install -m conda-forge cxx-compiler jupyter numpy matplotlib scipy sympy
```

Then launch Jupyter

```sh
jupyter notebook
```

and chose one of tutorial:

1. `lotka_volterra_rungekutta.ipynb`: first tutorial where explains how to start with Ponio with a Runge-Kutta method.
2. `lotka_volterra_lawson.ipynb`: how to use a Lawson method in Ponio.
3. `lotka_volterra_splitting.ipynb`: how to use a Lie or Strang splitting method.
4. `lotka_volterra_tuto.ipynb`: same as 3 firsts examples with tutorial of how to implement an observer
5. `arenstorf_adaptivetimestep.ipynb`: example with adaptive time step methods (and underlying Lawson methods)
