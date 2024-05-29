# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

name = "curtiss_hirschfelder_exprk"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

methods = {
    "rk44": "RK(4, 4)",
    "krogstad": "expRK(4, 4) Krogstad",
    "lrk44": "LRK(4, 4)"
}

plt.rcParams["figure.figsize"] = (8, 4)

line_style = "-"
for tag, meth in methods.items():
    data = np.loadtxt(os.path.join(data_dir, f"{tag}.dat"))
    t, y, dt = data[:, 0], data[:, 1], data[:, 2]

    plt.plot(t, y, line_style, label=meth)
    line_style = "--"

plt.ylabel("$y$")
plt.xlabel('time')

plt.legend()

plt.savefig(os.path.join(data_dir, "01.png"), dpi=200)

plt.title("Solution of Curtiss Hirschfelder problem")
plt.show()
