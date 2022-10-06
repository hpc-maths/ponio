#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
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

for i,x0 in enumerate(X):
  data = np.loadtxt(filename_fmt.format(x0))
  plt.plot(data[:,1],data[:,2],"-",color=cm.tab10(i),label=r"$x_0 = {:0.2f}$".format(x0))
  plt.plot([x0],[x0],"o",color=cm.tab10(i),label=None)

plt.title("Lotka-Volterra system")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.legend()
plt.show()
