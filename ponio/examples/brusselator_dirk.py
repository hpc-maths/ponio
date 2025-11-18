# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import dataclasses
import argparse

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

name = "brusselator_dirk"
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


@dataclasses.dataclass
class simu_data:
    label: str
    data: npt.NDArray[np.float64]

    def time(self):
        return self.data[:, 0]

    def concentration_x(self):
        return self.data[:, 1]

    def concentration_y(self):
        return self.data[:, 2]


def load_data(tagname, label):
    return simu_data(label=label, data=np.loadtxt(os.path.join(data_dir, f"brusselator_{tagname}.dat")))


def plot_concentration(data, line):
    plt.plot(data.time(), data.concentration_x(),
             line, label=f"x with {data.label}")
    plt.plot(data.time(), data.concentration_y(),
             line, label=f"y with {data.label}")


def plot_phase(data, line):
    plt.plot(data.concentration_x(),
             data.concentration_y(), line, label=data.label)


make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

methods = {
    'dirk23': "DIRK (2,3)",
    'dirk23_exact_solver': "DIRK (2,3) (exact solver)",
    'rk33': "RK (3,3)",
}

lines = ["+-", "x:", "+--"]

# concentration plot
for i, (tag, label) in enumerate(methods.items()):
    data = load_data(tag, label)
    plot_concentration(data, lines[i])

plt.xlabel('time')
plt.legend(loc="upper right")

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    plt.title("Solution of Brusselator system")
    plt.show()

# concentration plot in phase space
for i, (tag, label) in enumerate(methods.items()):
    data = load_data(tag, label)
    plot_phase(data, lines[i])

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc="upper right")

plt.savefig(os.path.join(img_dir, "02.png"), dpi=200)

if not arguments.only_save:
    plt.title("Phase plane")
    plt.show()
