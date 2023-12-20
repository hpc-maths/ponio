# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

name = "heat_rock"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt(os.path.join(data_dir, "heat_sol_rock2.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="solution with ROCK2")

data = np.loadtxt(os.path.join(data_dir, "heat_sol_rock4.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="solution with ROCK4")

data = np.loadtxt(os.path.join(data_dir, "heat_qexa.dat"))
plt.plot(data[:, 0], data[:, 1], ":",
         label="quasi-exact solution with RKC(20, 2)")

data = np.loadtxt(os.path.join(data_dir, "heat_ini.dat"))
plt.plot(data[:, 0], data[:, 1], ":", label="initial solution")

plt.title("Heat equation")
plt.legend(loc="lower left")
plt.ylim(bottom=-0.4)
plt.show()

data = np.loadtxt(os.path.join(data_dir, "errors.dat"))
plt.plot(data[:, 0], data[:, 1], "+-", label="ROCK2 errors in $\|\cdot\|_2$")
plt.plot(data[:, 0], data[:, 2], "+-", label="ROCK4 errors in $\|\cdot\|_2$")
plt.plot(data[:, 0], data[:, 0]**2, "--", label="slope order 2")
plt.plot(data[:, 0], (2.25*data[:, 0])**4, "--", label="slope order 4")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.show()
