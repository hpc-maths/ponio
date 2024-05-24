# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

name = "brusselator"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt(os.path.join(data_dir, "brusselator.dat"))
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

# concentration plot
plt.plot(t, x, "-", label="$x$")
plt.plot(t, y, "-", label="$y$")

plt.xlabel('time')
plt.legend()

plt.savefig("3-brusselator-equations_01.png", dpi=200)

plt.title("Solution of Brusselator system with RK (8,6) method")
plt.show()

# concentration plot in phase space
plt.plot(x, y, "-")
plt.xlabel('x')
plt.ylabel('y')

plt.savefig("3-brusselator-equations_02.png", dpi=200)

plt.title("Phase plane")
plt.show()
