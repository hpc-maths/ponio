#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
name = "arenstorf"
data_dir = name+"_data"

os.makedirs(data_dir, exist_ok=True)

make = subprocess.Popen(["make",name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data_1 = np.loadtxt(os.path.join(data_dir,"arenstorf_rk546m.dat"))
data_2 = np.loadtxt(os.path.join(data_dir,"arenstorf_rk547m.dat"))
data_3 = np.loadtxt(os.path.join(data_dir,"arenstorf_rk547s.dat"))

plt.plot(data_1[:,1],data_1[:,2],"-+",label="RK5(4) 6M")
plt.plot(data_2[:,1],data_2[:,2],"-+",label="RK5(4) 7M")
plt.plot(data_3[:,1],data_3[:,2],"-+",label="RK5(4) 7S")

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