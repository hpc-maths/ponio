# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import dataclasses
import glob

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import numpy.typing as npt

name = "lorenz_all"
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


def load_data(tagname):
    return simu_data(label=tagname, data=np.loadtxt(os.path.join(data_dir, f"{tagname}.dat")))


make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = [load_data(os.path.basename(f)[:-4])
        for f in glob.glob(os.path.join(data_dir, "*.dat"))]

# 3d plot
fig = plt.figure(constrained_layout=True, figsize=(8, 8))
gs = fig.add_gridspec(6, 1)
axlist = []

ax = fig.add_subplot(gs[:, 0], projection='3d')
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

lines = [ax.plot(d.x()[0], d.y()[0], d.z()[0], linewidth=0.375)[0]
         for d in data]


def update(frame):
    for l, d in zip(lines, data):
        l.set_data_3d(d.x()[:frame], d.y()[:frame], d.z()[:frame])

    return lines

# for i, d in enumerate(data):
#     ax.plot(d.x(), d.y(), d.z(), label=d.label,
#             color=f"#{((712*i) % (16**4)):04x}a2", linewidth=0.375)


ani = animation.FuncAnimation(
    fig=fig, func=update, frames=len(data[0].time()), interval=1)

# ax.legend()
axlist.append(ax)

print("save...")
ani.save(filename=os.path.join(data_dir, "01.gif"), writer="pillow")
print("!")
# plt.savefig(os.path.join(data_dir, "01.png"))

ax.set_title("Lorenz Attractor")
plt.show()
