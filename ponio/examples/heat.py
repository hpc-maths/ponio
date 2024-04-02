# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

name = "heat"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

plt.rcParams["figure.figsize"] = (8, 4)

data = np.loadtxt(os.path.join(data_dir, "heat_ini.dat"))
plt.plot(data[:, 0], data[:, 1], ":", label="initial condition")

data = np.loadtxt(os.path.join(data_dir, "heat_sol.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="solution with RKC (15,2)")

data = np.loadtxt(os.path.join(data_dir, "heat_exa.dat"))
plt.plot(data[:, 0], data[:, 1], ":", label="exact solution")

plt.legend()
plt.ylim(-0.1, 1)

plt.savefig(os.path.join(data_dir, "01.png"))

plt.title("Heat equation")
plt.show()
