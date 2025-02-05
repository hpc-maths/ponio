# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
import h5py


def import_2d_h5(filename_f5):
    mesh = h5py.File(filename_f5, 'r')['mesh']
    points = mesh['points']
    connectivity = mesh['connectivity']
    squares = np.array([points[connectivity[i]][:, :-1]
                       for i in range(connectivity.shape[0])])
    centers = np.array([0.25*np.sum(pts, axis=0) for pts in squares])

    return squares, centers, mesh['fields']


class bz_data:
    def __init__(self, filename_h5):
        self.squares, self.centers, self.fields = import_2d_h5(filename_h5)

    def x(self):
        return self.centers[:, 0][:]

    def y(self):
        return self.centers[:, 1][:]

    def z(self, field):
        return self.fields[field][:]

    def dx(self):
        return np.max(self.squares[:, :, 0], axis=1) - np.min(self.squares[:, :, 0], axis=1)

    def dy(self):
        return np.max(self.squares[:, :, 1], axis=1) - np.min(self.squares[:, :, 1], axis=1)

    def dxdy(self):
        return self.dx()*self.dy()

    def integral(self, field):
        return np.sum(self.z(field)*self.dxdy())


name = "bz_2d_pirock"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

print(os.path.join(data_dir, "y_final.h5"))
data = bz_data(os.path.join(data_dir, "y_final.h5"))

# b and c concentration plot
fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.tricontour(data.x(), data.y(), data.z(
    "u_1"), levels=5, linewidths=0.5, colors='k')
ax.tricontourf(data.x(), data.y(), data.z(
    "u_2"), levels=200, cmap="RdBu_r")

ax.plot(np.inf, np.inf, "-", linewidth=0.5, color='k', label="$b$")
ax.plot(np.inf, np.inf, "s", color='red', alpha=0.75, label="$c$")
ax.legend()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")

plt.savefig(os.path.join(data_dir, "01.png"), dpi=200)

plt.show()

# levels plot
fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.tricontour(data.x(), data.y(), data.z("level"), levels=range(
    min(data.z("level")), max(data.z("level"))+1), linewidths=0.5, colors='k')
ax.tricontourf(data.x(), data.y(), data.z(
    "u_0"), levels=200, cmap="RdBu_r")

ax.plot(np.inf, np.inf, "-", linewidth=0.5, color='k', label="levels")
ax.plot(data.x(), data.y(), "+", color="k",
        alpha=0.125, label="center of cells")
ax.legend()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")

plt.savefig(os.path.join(data_dir, "02.png"), dpi=200)

plt.show()
