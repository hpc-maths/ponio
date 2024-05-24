# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

name = "pendulum"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt(os.path.join(data_dir, "pendulum.dat"))


def pendulum(t, y, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt


y0 = [np.pi-0.1, 0.0]
t_span = (0, 10)
t = np.linspace(0, 10, 101)
b = 0.25
c = 5.0
sol = solve_ivp(pendulum, t_span, y0, t_eval=t, args=(b, c))

plt.plot(t, sol.y[0], '-x', label=r"$\theta$ with SciPy (solve_ivp)")
plt.plot(t, sol.y[1], '-x', label=r"$\omega$ with SciPy (solve_ivp)")

plt.plot(data[:, 0], data[:, 1], "-+", label=r"$\theta$ with ponio")
plt.plot(data[:, 0], data[:, 2], "-+", label=r"$\omega$ with ponio")

plt.xlabel("time")
plt.legend()

plt.savefig("15-pendulum-equation_01.png")

plt.title("Pendulum equation (solved with RK (4,4))")
plt.show()
