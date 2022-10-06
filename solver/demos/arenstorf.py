#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

os.makedirs("example4",exist_ok=True)

make = subprocess.Popen(["make","arenstorf"])
make.wait()

args = ["./arenstorf"]
process = subprocess.Popen(args)
process.wait()

data_1 = np.loadtxt("example4/arenstorf_rk546m.dat")
data_2 = np.loadtxt("example4/arenstorf_rk547m.dat")
data_3 = np.loadtxt("example4/arenstorf_rk547s.dat")

plt.plot(data_1[:,1],data_1[:,3],"-+",label="RK5(4) 6M")
plt.plot(data_2[:,1],data_2[:,3],"-+",label="RK5(4) 7M")
plt.plot(data_3[:,1],data_3[:,3],"-+",label="RK5(4) 7S")

plt.axhline(0.,-1.3,1.3,linestyle="--",color="grey",linewidth=1)
plt.axvline(0.,-1.1,1.1,linestyle="--",color="grey",linewidth=1)

plt.axis('equal')

plt.title("Arenstorf orbit")
plt.legend()
plt.show()

plt.title("Arenstorf orbit, time step size")
plt.plot(data_1[:,0],data_1[:,-1],"-+",label="RK5(4) 6M")
plt.plot(data_2[:,0],data_2[:,-1],"-+",label="RK5(4) 7M")
plt.plot(data_3[:,0],data_3[:,-1],"-+",label="RK5(4) 7S")

plt.xlabel(r"$t$")
plt.ylabel(r"$\Delta t$")
plt.legend()
plt.show()