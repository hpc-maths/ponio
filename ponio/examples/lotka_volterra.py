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

name = "lotka_volterra"
data_dir = name+"_data"
filename_fmt = f"{name}_{{:0.2f}}.dat"
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
    x0: float
    data: npt.NDArray[np.float64]

    def time(self):
        return self.data[:, 0]

    def rabbit(self):
        return self.data[:, 1]

    def fox(self):
        return self.data[:, 2]


def load_data(x0):
    return simu_data(x0=x0, data=np.loadtxt(os.path.join(data_dir, filename_fmt.format(x0))))


os.makedirs(data_dir, exist_ok=True)

make = subprocess.Popen(["make", name])
make.wait()

X0 = np.linspace(0.9, 1.8, 10)

# launch multiple simulation
for x0 in X0:
    args = [os.path.join(".", name),
            os.path.join(data_dir, filename_fmt.format(x0)), str(x0)]
    process = subprocess.Popen(args)
    process.wait()

data = [load_data(x0) for x0 in X0]

# composantes plot
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 6))

for i, d in enumerate(data):
    ax1.plot(d.time(), d.rabbit(), color=f"C{i}", label=f"$x_0 = {d.x0:0.2f}$")
    ax2.plot(d.time(), d.fox(), color=f"C{i}", label=f"$x_0 = {d.x0:0.2f}$")

ax1.set_ylabel("prey: $x$")
ax2.set_ylabel("predator: $y$")
ax2.set_xlabel("time")

box1 = ax1.get_position()
ax1.set_position([box1.x0, box1.y0, box1.width * 0.95, box1.height])
box2 = ax2.get_position()
ax2.set_position([box2.x0, box2.y0, box2.width * 0.95, box2.height])

ax1.legend(bbox_to_anchor=(1.01, 1),
           loc='upper left', borderaxespad=0.)

plt.savefig(os.path.join(img_dir, "01.png"))

if not arguments.only_save:
    fig.suptitle("Lotka-Volterra system (solved witk RK (11,8))")
    plt.show()


for i, d in enumerate(data):
    plt.plot(d.rabbit(), d.fox(), color=f"C{i}", label=f"$x_0 = {d.x0:0.2f}$")
    plt.plot([d.x0], [d.x0], "o", color=f"C{i}", label=None)

plt.xlabel("$x$")
plt.ylabel("$y$")
plt.legend()

plt.savefig(os.path.join(img_dir, "02.png"))

if not arguments.only_save:
    plt.title("Lotka-Volterra system (solved witk RK (11,8))")
    plt.show()
