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

#####################################################################
# low def 3d plot
fig = plt.figure(figsize=(6, 4))
fig.subplots_adjust(top=1.25, bottom=-0.25, left=-0.5, right=1.5)
gs = fig.add_gridspec(1, 1)
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
plt.axis('off')
ax.grid(False)

lines = [ax.plot(d.x()[0], d.y()[0], d.z()[0], linewidth=0.375)[0]
         for d in data]

N_frames = 42


def update(frame):
    for l, d in zip(lines, data):
        len_simu = len(d.time())
        i = int(frame/N_frames*len_simu)
        l.set_data_3d(d.x()[:i], d.y()[:i], d.z()[:i])

    return lines


ani = animation.FuncAnimation(
    fig=fig, func=update, frames=N_frames, interval=1)

axlist.append(ax)

print("save...")
ani.save(filename=os.path.join(data_dir, "01.gif"), dpi=50,
         writer="pillow")
print("!")

# ax.set_title("Lorenz Attractor")
# plt.show()

#####################################################################
# full 3d plot
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(top=1., bottom=0., left=-0.5, right=1.5)
gs = fig.add_gridspec(1, 1)
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
plt.axis('off')
ax.grid(False)

lines = [ax.plot(d.x()[0], d.y()[0], d.z()[0], linewidth=0.375, color=f"C{i}")[0]
         for i, d in enumerate(data)]
points = [ax.plot(d.x()[0], d.y()[0], d.z()[0], "o", linewidth=0.375, color=f"C{i}")[0]
          for i, d in enumerate(data)]

N_frames = 256


def update(frame):
    for p, l, d in zip(points, lines, data):
        len_simu = len(d.time())
        i = int(frame/N_frames*len_simu)
        j = 0  # max(0, i-200)
        p.set_data_3d([d.x()[i]], [d.y()[i]], [d.z()[i]])
        l.set_data_3d(d.x()[j:i], d.y()[j:i], d.z()[j:i])

    return lines


ani = animation.FuncAnimation(
    fig=fig, func=update, frames=N_frames, interval=1)

axlist.append(ax)

ax.set_title("Lorenz Attractor")
plt.show()
