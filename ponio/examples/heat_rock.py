# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

name = "heat_rock"
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

data = np.loadtxt(os.path.join(data_dir, "heat_sol_rock2.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="solution with ROCK2")

data = np.loadtxt(os.path.join(data_dir, "heat_sol_rock4.dat"))
plt.plot(data[:, 0], data[:, 1], "-", label="solution with ROCK4")

data = np.loadtxt(os.path.join(data_dir, "heat_qexa.dat"))
plt.plot(data[:, 0], data[:, 1], ":",
         label="quasi-exact solution with RKC(20, 2)")

data = np.loadtxt(os.path.join(data_dir, "heat_ini.dat"))
plt.plot(data[:, 0], data[:, 1], ":", label="initial solution")

plt.legend(loc="lower left")
plt.ylim(bottom=-0.4)

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    plt.title("Heat equation")
    plt.show()

data = np.loadtxt(os.path.join(data_dir, "errors.dat"))
plt.plot(data[:, 0], data[:, 1], "+-", label="ROCK2 errors in $\|\cdot\|_2$")
plt.plot(data[:, 0], data[:, 2], "+-", label="ROCK4 errors in $\|\cdot\|_2$")
plt.plot(data[:, 0], data[:, 0]**2, "--", label="slope order 2")
plt.plot(data[:, 0], (2.25*data[:, 0])**4, "--", label="slope order 4")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()

plt.savefig(os.path.join(img_dir, "02.png"), dpi=200)

if not arguments.only_save:
    plt.show()
