# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import dataclasses

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

name = "arenstorf"
data_dir = name+"_data"


@dataclasses.dataclass
class simu_data:
    label: str
    data: npt.NDArray[np.float64]

    def time(self):
        return self.data[:, 0]

    def pos_x(self):
        return self.data[:, 1]

    def pos_y(self):
        return self.data[:, 2]

    def vel_x(self):
        return self.data[:, 3]

    def vel_y(self):
        return self.data[:, 4]

    def time_step(self):
        return self.data[:, -1]


def load_data(tagname, label):
    return simu_data(label=label, data=np.loadtxt(os.path.join(data_dir, f"{name}_{tagname}.dat")))


def plot_orbit(data):
    plt.plot(data.pos_x(), data.pos_y(), "-+", label=data.label)


def plot_velocity(data):
    plt.plot(data.vel_x(), data.vel_y(), "-+", label=data.label)


def plot_time_step(data):
    plt.plot(data.time(), data.time_step(), "-+", label=data.label)


os.makedirs(data_dir, exist_ok=True)

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

rk546m = load_data("rk546m", "RK5(4) 6M")
rk547m = load_data("rk547m", "RK5(4) 7M")
rk547s = load_data("rk547s", "RK5(4) 7S")

# orbit plot
plt.rcParams["figure.figsize"] = (6, 6)
plt.axes().set_aspect('equal', 'box')

plot_orbit(rk546m)
plot_orbit(rk547m)
plot_orbit(rk547s)

plt.axhline(0., 0, 1, linestyle="--", color="grey", linewidth=1)
plt.axvline(0., 0, 1, linestyle="--", color="grey", linewidth=1)

plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)

plt.xlabel("$x$")
plt.ylabel("$y$")

plt.tight_layout(pad=1.0)
plt.legend()

plt.savefig("1-arenstorf-orbit_01.png", dpi=200)

plt.title("Arenstorf orbit")
plt.show()

plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]

# velocity plot
plt.rcParams["figure.figsize"] = (6, 6)
plt.axes().set_aspect('equal', 'box')

plot_velocity(rk546m)
plot_velocity(rk547m)
plot_velocity(rk547s)

plt.axhline(0., 0, 1, linestyle="--", color="grey", linewidth=1)
plt.axvline(0., 0, 1, linestyle="--", color="grey", linewidth=1)

plt.xlim(-2.5, 2.5)
plt.ylim(-3.5, 1.5)

plt.xlabel(r"$\dot{{x}}$")
plt.ylabel(r"$\dot{{y}}$")

plt.tight_layout(pad=1.0)
plt.legend()

plt.savefig("1-arenstorf-orbit_02.png", dpi=200)

plt.title("Arenstorf velocity")
plt.show()

plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]

# time step history
plot_time_step(rk546m)
plot_time_step(rk547m)
plot_time_step(rk547s)

plt.xlabel(r"$t$")
plt.ylabel(r"$\Delta t$")
plt.legend()

plt.savefig("1-arenstorf-orbit_03.png")

plt.title("Arenstorf orbit, time step size")
plt.show()
