# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
import h5py


def import_h5(filename_f5):
    file = h5py.File(filename_f5, 'r')
    mesh = file['mesh']
    points = mesh['points']
    connectivity = mesh['connectivity']
    squares = np.array([points[connectivity[i]][:, :-1]
                       for i in range(connectivity.shape[0])])
    centers = np.array([np.mean(pts, axis=0) for pts in squares])
    fields = {
        key: mesh['fields'][key][:]
        for key in mesh['fields'].keys()
    }
    file.close()

    return squares, centers, fields


class h5_data:
    def __init__(self, filename_h5):
        self.squares, self.centers, self.fields = import_h5(filename_h5)

    @property
    def x(self):
        return self.centers[:, 0][:]

    @property
    def y(self):
        return self.centers[:, 1][:]

    @property
    def z(self):
        return self.centers[:, 1][:]

    def __getitem__(self, field):
        return self.fields[field][:]

    def keys(self):
        return self.fields.keys()


name = "brusselator_advection_2d"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = h5_data(os.path.join(data_dir, "uv_final.h5"))

# b and c concentration plot
fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.tricontour(data.x, data.y, data["u_0"],
              levels=5, linewidths=0.5, colors='k')
ax.tricontourf(data.x, data.y, data["u_1"], levels=200, cmap="RdBu_r")

ax.plot(np.inf, np.inf, "s", color='red', alpha=0.75, label="$b$")
ax.plot(np.inf, np.inf, "-", linewidth=0.5, color='k', label="$c$")
ax.legend()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")

plt.savefig(os.path.join(data_dir, "01.png"), dpi=200)

plt.show()

# levels plot
fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.tricontour(data.x, data.y, data["level"], levels=range(
    min(data["level"]), max(data["level"])+1), linewidths=0.5, colors='k')
ax.tricontourf(data.x, data.y, data["u_0"], levels=200, cmap="RdBu_r")

ax.plot(np.inf, np.inf, "-", linewidth=0.5, color='k', label="levels")
ax.plot(data.x, data.y, "+", color="k",
        alpha=0.125, label="center of cells")
ax.legend()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")

plt.savefig(os.path.join(data_dir, "02.png"), dpi=200)

plt.show()
