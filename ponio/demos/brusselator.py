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
t  = data[:,0]
u1 = data[:,1]
u2 = data[:,2]
plt.plot(t, u1,"-+", label="$u_1$")
plt.plot(t, u2,"-+", label="$u_2$")

plt.title("Solution of Brusselator system with RK (8,6) method")
plt.xlabel('time')
plt.legend()
plt.show()

plt.title("Phase plane")
plt.plot(u1, u2,"-+")
plt.xlabel('u1')
plt.ylabel('u2')
plt.show()
