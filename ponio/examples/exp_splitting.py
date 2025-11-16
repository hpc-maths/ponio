# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

name = "exp_splitting"
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

test_cases = [
    ("strang", "Strang splitting with exact solver and RK(2,2) Ralston"),
    ("exact", "exact solver"),
    ("rk2", "RK(2,2) Ralston")
]

for tag, label in test_cases:
    data = np.loadtxt(os.path.join(data_dir, f"exp_{tag}.dat"))

    plt.plot(data[:, 0], data[:, 1], "-+", label=label)

plt.plot(data[:, 0], np.exp(data[:, 0]), "--", label="exact solution")
plt.legend()

plt.savefig(os.path.join(img_dir, "01.png"))

if not arguments.only_save:
    plt.title(r"$\dot{u} = \lambda u + (1-\lambda) u$, $u(t=0)=1$")
    plt.show()
