#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

X = np.linspace(0.9,1.8,10)
filename_fmt = "lv/lokta_voltera_{:0.2f}.dat"

os.makedirs("lv",exist_ok=True)

subprocess.Popen(["make","lokta_voltera"])

for x0 in X:
  args = ["./lokta_voltera",filename_fmt.format(x0),str(x0)]
  subprocess.Popen(args)

for x0 in X:
  data = np.loadtxt(filename_fmt.format(x0))
  plt.plot(data[:,1],data[:,2],"-",label=r"$x_0 = {:0.2f}$".format(x0))

plt.title("Lokta-Voltera system")
plt.legend()
plt.show()
