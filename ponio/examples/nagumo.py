# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import dataclasses
import glob
import argparse

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

name = "nagumo"
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

k = 1.0
d = 1.0

x_0 = -10.


@dataclasses.dataclass
class simu_data:
    time: float
    data: npt.NDArray[np.float64]

    def x(self):
        return self.data[:, 0]

    def u(self):
        return self.data[:, 1]

    def exact_sol(self):
        v = (1./np.sqrt(2.))*np.sqrt(k*d)
        cst = -(1./np.sqrt(2.))*np.sqrt(k/d)

        return np.exp(cst*(self.x() - x_0 - v*self.time)) / (1. + np.exp(cst*(self.x() - x_0 - v*self.time)))


def load_data(filename):
    time = float(filename.split("/")[-1][:-4].split("_")[1])
    return simu_data(time=time, data=np.loadtxt(filename))


os.makedirs(data_dir, exist_ok=True)

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = sorted([load_data(file) for file in glob.glob(
    os.path.join(data_dir, "*.dat"))], key=lambda d: d.time)

for i, dat in enumerate(data):
    plt.plot(dat.x(), dat.u(), color=f"C{i}", label=f"$t^n = {dat.time}$")

plt.legend()

plt.savefig(os.path.join(img_dir, "01.png"))

if not arguments.only_save:
    plt.title("Nagumo solution")
    plt.show()

for i, dat in enumerate(data):
    plt.plot(dat.x(), dat.u() - dat.exact_sol(),
             color=f"C{i}", label=f"$t^n = {dat.time}$")

plt.legend()

plt.savefig(os.path.join(img_dir, "02.png"))

if not arguments.only_save:
    plt.title("Absolute error in time")
    plt.show()
