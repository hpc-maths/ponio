# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np

name = "brusselator_dirk"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

methods = {
            'dirk23': "DIRK (2,3)",
            'rk33': "RK (3,3)"
        }

for tag, label in methods.items():
    data = np.loadtxt(os.path.join(data_dir, f"brusselator_{tag}.dat"))
    t  = data[:,0]
    u1 = data[:,1]
    u2 = data[:,2]
    plt.plot(t, u1,"-+", label=f"$u_1$ with {label}")
    plt.plot(t, u2,"-+", label=f"$u_2$ with {label}")

plt.title("Solution of Brusselator system")
plt.xlabel('time')
plt.legend()
plt.show()

for tag, label in methods.items():
    data = np.loadtxt(os.path.join(data_dir, f"brusselator_{tag}.dat"))
    t  = data[:,0]
    u1 = data[:,1]
    u2 = data[:,2]
    plt.plot(u1, u2, "-+", label=label)

plt.title("Phase plane")
plt.xlabel('$u_1$')
plt.ylabel('$u_2$')
plt.legend()
plt.show()
