#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

os.makedirs("example3",exist_ok=True)

make = subprocess.Popen(["make","pendulum"])
make.wait()

args = ["./pendulum"]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt("example3/pendulum.dat")

def pendulum(y,t,b,c):
  theta,omega = y
  dydt = [omega, -b*omega - c*np.sin(theta)]
  return dydt
y0 = [np.pi-0.1,0.0]
t = np.linspace(0,10,101)
b = 0.25
c = 5.0
sol = odeint(pendulum,y0,t,args=(b,c))

plt.plot(t, sol[:, 0], '-x', label=r"$\theta$ with SciPy odeint")
plt.plot(t, sol[:, 1], '-x', label=r"$\theta$ with SciPy odeint")

plt.plot(data[:,0],data[:,1],"-+",label=r"$\theta$ with ponio")
plt.plot(data[:,0],data[:,2],"-+",label=r"$\omega$ with ponio")

plt.title("Pendulum equation")
plt.legend()
plt.show()
