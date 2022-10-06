#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

X = np.linspace(0.9,1.8,10)
filename_fmt = "example2/lotka_volterra_{:0.2f}.dat"

os.makedirs("example2",exist_ok=True)

make = subprocess.Popen(["make","lotka_volterra"])
make.wait()

for x0 in X:
  args = ["./lotka_volterra",filename_fmt.format(x0),str(x0)]
  process = subprocess.Popen(args)
process.wait()

for x0 in X:
  data = np.loadtxt(filename_fmt.format(x0))
  plt.plot(data[:,1],data[:,2],"-",label=r"$x_0 = {:0.2f}$".format(x0))

plt.title("Lotka-Volterra system")
plt.legend()
plt.show()
