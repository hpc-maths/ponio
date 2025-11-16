# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

name = "brusselator"
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

data = np.loadtxt(os.path.join(data_dir, "brusselator.dat"))
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

# concentration plot
plt.plot(t, x, "-", label="$x$")
plt.plot(t, y, "-", label="$y$")

plt.xlabel('time')
plt.legend()

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    plt.title("Solution of Brusselator system with RK (8,6) method")
    plt.show()

# concentration plot in phase space
plt.plot(x, y, "-")
plt.xlabel('x')
plt.ylabel('y')

plt.savefig(os.path.join(img_dir, "02.png"), dpi=200)

if not arguments.only_save:
    plt.title("Phase plane")
    plt.show()
