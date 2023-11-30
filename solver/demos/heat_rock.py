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

data = np.loadtxt(os.path.join(data_dir, "heat_sol.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="solution with ROCK")

data = np.loadtxt(os.path.join(data_dir, "heat_qexa.dat"))
plt.plot(data[:, 0], data[:, 1], ":",
         label="quasi-exact solution with RKC(20, 2)")

data = np.loadtxt(os.path.join(data_dir, "heat_ini.dat"))
plt.plot(data[:, 0], data[:, 1], ":", label="initial solution")

plt.title("Heat equation")
plt.legend()
# plt.ylim(top=0.8)
plt.show()

data = np.loadtxt(os.path.join(data_dir, "errors.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="errors in $\|\cdot\|_2$")
plt.plot(data[:, 0], data[:, 0]**2, "--", label="slope order 2")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.show()
