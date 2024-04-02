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

name = "lorenz_tuto"
data_dir = name+"_data"


@dataclasses.dataclass
class simu_data:
    label: str
    data: npt.NDArray[np.float64]

    def time(self):
        return self.data[:, 0]

    def x(self):
        return self.data[:, 1]

    def y(self):
        return self.data[:, 2]

    def z(self):
        return self.data[:, 3]


def load_data(tagname, label):
    return simu_data(label=label, data=np.loadtxt(os.path.join(data_dir, f"{tagname}.dat")))


make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()


rk44 = load_data("rk44", "RK(4, 4)")
lrk44 = load_data("lrk44", "LRK(4, 4)")
lie = load_data("lie", "Lie")
strang = load_data("strang", "Strang")

data = [load_data(tag, label) for (tag, label) in [
    ("rk44", "RK(4,4)"), ("lrk44", "LRK(4,4)"), ("lie", "Lie"), ("strang", "Strang")]]

# 3d plot
fig = plt.figure(constrained_layout=True, figsize=(8, 8))
gs = fig.add_gridspec(6, 1)
axlist = []

ax = fig.add_subplot(gs[:, 0], projection='3d')
ax.set_title("Lorenz attractor")
ax.set_xlim3d([-20, 20])
ax.set_ylim3d([-20, 20])
ax.set_zlim3d(bottom=0, top=50)

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')
# ax.grid(False)

for d in data:
    ax.plot(d.x(), d.y(), d.z(), label=d.label, linewidth=0.375)

ax.legend()
axlist.append(ax)

plt.savefig(os.path.join(data_dir, "01.png"))

ax.set_title("Lorenz Attractor")
plt.show()

# 2d plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)

for d in data:
    ax1.plot(d.time(), d.x(), linewidth=0.5, label=d.label)

ax1.set_ylabel("$x$")
xticklabels = ax1.get_xticklabels()
for xlab in xticklabels:
    xlab.set_text("")
ax1.set_xticks(ax1.get_xticks(), xticklabels)
ax1.set_xlim(d.time()[0], d.time()[-1])
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncols=4, mode="expand", borderaxespad=0.)

for d in data:
    ax2.plot(d.time(), d.y(), linewidth=0.5, label=d.label)

ax2.set_ylabel("$y$")
ax2.set_xticks(ax2.get_xticks(), xticklabels)
ax2.set_xlim(d.time()[0], d.time()[-1])

for d in data:
    ax3.plot(d.time(), d.z(), linewidth=0.5, label=d.label)

ax3.set_ylabel("$z$")

ax3.set_xlabel("$t$")
ax3.set_xlim(d.time()[0], d.time()[-1])

plt.savefig(os.path.join(data_dir, "01.png"))

fig.suptitle("Composantes of Lorenz Attractor")
plt.show()
