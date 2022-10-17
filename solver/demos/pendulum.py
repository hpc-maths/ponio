#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

os.makedirs("example3",exist_ok=True)

make = subprocess.Popen(["make","pendulum"])
make.wait()

args = ["./pendulum"]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt("example3/pendulum.dat")

def pendulum(t,y,b,c):
  theta,omega = y
  dydt = [omega, -b*omega - c*np.sin(theta)]
  return dydt
y0 = [np.pi-0.1,0.0]
t_span = (0,10)
t = np.linspace(0,10,101)
b = 0.25
c = 5.0
sol = solve_ivp(pendulum,t_span,y0,t_eval=t,args=(b,c))

plt.plot(t, sol.y[0], '-x', label=r"$\theta$ with SciPy (solve_ivp)")
plt.plot(t, sol.y[1], '-x', label=r"$\omega$ with SciPy (solve_ivp)")

plt.plot(data[:,0],data[:,1],"-+",label=r"$\theta$ with ponio")
plt.plot(data[:,0],data[:,2],"-+",label=r"$\omega$ with ponio")

plt.title("Pendulum equation (solved with RK (4,4))")
plt.legend()
plt.show()