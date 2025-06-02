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

    def norm_l1(self):
        return np.abs(self.data[:, 1]) + np.abs(self.data[:, 2]) + np.abs(self.data[:, 3]) + np.abs(self.data[:, 4])

    def norm_l2(self):
        return np.sqrt(np.square(self.data[:, 1]) + np.square(self.data[:, 2]) + np.square(self.data[:, 3]) + np.square(self.data[:, 4]))


def load_data(tagname, label):
    return simu_data(label=label, data=np.loadtxt(os.path.join(data_dir, f"{name}_{tagname}.dat")))


def plot_orbit(data, alpha=1):
    plt.plot(data.pos_x(), data.pos_y(), "-+", alpha=alpha, label=data.label)


def plot_velocity(data, alpha=1):
    plt.plot(data.vel_x(), data.vel_y(), "-+", alpha=alpha, label=data.label)


def plot_time_step(data, alpha=1):
    plt.plot(data.time(), data.time_step(),
             "-+", alpha=alpha, label=data.label)


def plot_error(data, ref, alpha=1):
    plt.plot(data.time(), np.abs(data.norm_l2() - np.interp(data.time(),
             ref.time(), ref.norm_l2())), "-+", alpha=alpha, label=data.label)


os.makedirs(data_dir, exist_ok=True)

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

rk118 = load_data("rk118", "RK(11, 8)")
rk546m = load_data("rk546m", "RK5(4) 6M")
rk547m = load_data("rk547m", "RK5(4) 7M")
rk547s = load_data("rk547s", "RK5(4) 7S")
strang = load_data("adaptive_strang", "Adaptive Strang")

# orbit plot
plt.rcParams["figure.figsize"] = (6, 6)
plt.axes().set_aspect('equal', 'box')

plot_orbit(rk546m)
plot_orbit(rk547m)
plot_orbit(rk547s)
plot_orbit(strang, 0.05)

plt.axhline(0., 0, 1, linestyle="--", color="grey", linewidth=1)
plt.axvline(0., 0, 1, linestyle="--", color="grey", linewidth=1)

plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)

plt.xlabel("$x$")
plt.ylabel("$y$")

plt.tight_layout(pad=1.0)
plt.legend()

plt.savefig(os.path.join(data_dir, "01.png"), dpi=200)

plt.title("Arenstorf orbit")
plt.show()

plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]

# velocity plot
plt.rcParams["figure.figsize"] = (6, 6)
plt.axes().set_aspect('equal', 'box')

plot_velocity(rk546m)
plot_velocity(rk547m)
plot_velocity(rk547s)
plot_velocity(strang, 0.05)

plt.axhline(0., 0, 1, linestyle="--", color="grey", linewidth=1)
plt.axvline(0., 0, 1, linestyle="--", color="grey", linewidth=1)

plt.xlim(-2.5, 2.5)
plt.ylim(-3.5, 1.5)

plt.xlabel(r"$\dot{{x}}$")
plt.ylabel(r"$\dot{{y}}$")

plt.tight_layout(pad=1.0)
plt.legend()

plt.savefig(os.path.join(data_dir, "02.png"), dpi=200)

plt.title("Arenstorf velocity")
plt.show()

plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]

# time step history
plot_time_step(rk546m)
plot_time_step(rk547m)
plot_time_step(rk547s)
plot_time_step(strang)

plt.yscale('log')
plt.xlabel(r"$t$")
plt.ylabel(r"$\Delta t$")
plt.legend()

plt.savefig(os.path.join(data_dir, "03.png"))

plt.title("Arenstorf orbit, time step size")
plt.show()

# error throw time
plot_error(rk546m, rk118)
plot_error(rk547m, rk118)
plot_error(rk547s, rk118)
plot_error(strang, rk118, 0.1)

plt.yscale('log')
plt.xlabel(r"$t$")
plt.ylabel(r"$\|\cdot\|_2$")
plt.legend()

plt.savefig(os.path.join(data_dir, "04.png"))

plt.title("Arenstorf orbit, error over time")
plt.show()
