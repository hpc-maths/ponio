# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

name = "exp"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt(os.path.join(data_dir, "exp.dat"))

plt.plot(data[:, 0], np.exp(data[:, 0]), "--", label="exact solution")
plt.plot(data[:, 0], data[:, 1], "-+", label="solution with RK NSSP (2,1)")

plt.legend()

plt.savefig(os.path.join(data_dir, "01.png"))

plt.title("$\dot{u} = u$, $u(t=0)=1$")
plt.show()
