#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

os.makedirs("example1",exist_ok=True)

make = subprocess.Popen(["make","exp"])
make.wait()

args = ["./exp"]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt("example1/exp.dat")
plt.plot(data[:,0],np.exp(data[:,0]),"--",label="exact solution")
plt.plot(data[:,0],data[:,1],"-+",label="solution with RK NSSP (2,1)")

plt.title("$\dot{u} = u$, $u(t=0)=1$")
plt.legend()
plt.show()
