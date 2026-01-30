# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import argparse

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

name = "curtiss_hirschfelder"
data_dir = "ch_data"
img_dir = data_dir

parser = argparse.ArgumentParser(
    description=f"Compile, launch and plot results of example `{name}`")
parser.add_argument('--only-save', action='store_true',
                    help="Just save output in img directory")

arguments = parser.parse_args()

if arguments.only_save:
    img_dir = os.path.join("img", name)
    os.makedirs(img_dir, exist_ok=True)


make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

methods = {
    "rk_33_ralston_cst": "RK(3,3) Ralston fixed $\\Delta t$",
    "dirk23_cst": "DIRK(2, 3) fixed $\\Delta t$",
    "rk_33_ralston": "RK(3,3) Ralston",
    "rk54_6m": "RK5(4) 6m"
}

fig, axs = plt.subplots(2)
# axs[0].axvline(0.464, linewidth=1, linestyle="--", color="grey")
# axs[1].axvline(0.464, linewidth=1, linestyle="--", color="grey")

line_style = "+-"
for tag, meth in methods.items():
    data = np.loadtxt(os.path.join(data_dir, f"sol_{tag}.dat"))
    t, y, dt = data[:, 0], data[:, 1], data[:, 2]

    axs[0].plot(t, y, line_style, label=f"{meth}")
    axs[1].plot(t, dt, line_style, label=f"{meth}")
    line_style = "x-"

axs[0].set_ylabel("$y$")
axs[1].set_ylabel("$\\Delta t$")
axs[1].set_xlabel('time')

plt.legend(loc="lower right")

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    fig.suptitle("Solution of Curtiss Hirschfelder problem")
    plt.show()
