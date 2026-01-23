# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

name = "exp"
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

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt(os.path.join(data_dir, "exp.dat"))

plt.plot(data[:, 0], np.exp(data[:, 0]), "--", label="exact solution")
plt.plot(data[:, 0], data[:, 1], "-+", label="solution with RK NSSP (2,1)")

plt.legend()

plt.savefig(os.path.join(img_dir, "01.png"))

if not arguments.only_save:
    plt.title("$\\dot{u} = u$, $u(t=0)=1$")
    plt.show()
