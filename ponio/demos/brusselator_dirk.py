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

name = "brusselator_dirk"
data_dir = name+"_data"


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

plt.savefig("4-brusselator-equations-with-dirk-method_01.png", dpi=200)

plt.title("Solution of Brusselator system")
plt.show()

# concentration plot in phase space
for i, (tag, label) in enumerate(methods.items()):
    data = load_data(tag, label)
    plot_phase(data, lines[i])

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc="upper right")

plt.savefig("4-brusselator-equations-with-dirk-method_02.png", dpi=200)

plt.title("Phase plane")
plt.show()
