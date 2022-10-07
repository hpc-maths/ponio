#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

os.makedirs("example5",exist_ok=True)

make = subprocess.Popen(["make","brownian"])
make.wait()

N = 10

args = ["./brownian",str(N)]
process = subprocess.Popen(args)
process.wait()

for i in range(N):
  data = np.loadtxt("example5/brownian_{}.dat".format(i))
  plt.plot(data[:,1],data[:,2],"-",linewidth=1)

plt.axis('equal')

plt.title("Brownian motion")
plt.show()
