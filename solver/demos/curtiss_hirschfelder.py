# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

name = "curtiss_hirschfelder"
data_dir = "ch_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

for meth in ("rk_33_ralston", "rk54_6m"):
    data = np.loadtxt(os.path.join(data_dir, f"sol_{meth}.dat"))
    t, y = data[:,0], data[:,1]
    plt.plot(t, y, "-+", label=f"{meth} : $t_{{final}} = {t[-1]}$")

plt.title("Solution of Curtiss Hirschfelder problem")
plt.xlabel('time')
plt.legend()
plt.show()
