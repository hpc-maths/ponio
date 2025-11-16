# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

name = "brownian"
data_dir = name+"_data"
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

n = 10

args = [os.path.join(".", name), str(n)]
process = subprocess.Popen(args)
process.wait()

plt.rcParams["figure.figsize"] = (6, 6)
plt.axes().set_aspect('equal', 'box')

for i in range(n):
    data = np.loadtxt(os.path.join(data_dir, "brownian_{}.dat".format(i)))
    plt.plot(data[:, 1], data[:, 2], "-", linewidth=1)

plt.xlim(-.5, .5)
plt.ylim(-.5, .5)

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    plt.title("Brownian motion")
    plt.show()
