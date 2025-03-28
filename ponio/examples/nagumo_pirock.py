# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import glob

import matplotlib.pyplot as plt
import numpy as np
import h5py

name = "nagumo_pirock"
data_dir = name+"_data"


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

    def __getitem__(self, field):
        return self.fields[field][:]

    def keys(self):
        return self.fields.keys()


def exact_sol_gen(d, k, x_0):
    v = (1./np.sqrt(2.))*np.sqrt(k*d)
    cst = -(1./np.sqrt(2.))*np.sqrt(k/d)

    return lambda t, x: np.exp(cst*(x - x_0 - v*t)) / (1. + np.exp(cst*(x - x_0 - v*t)))


d = .1
k = 1./d

x_0 = -25.


os.makedirs(data_dir, exist_ok=True)

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

# assume we want read all h5 files in this directory and they all finish with `***_n.h5` with `n` the number of iteration
# 1. list all "*.h5" files in given directory
# 2. filter data which doesn't finish with `_n.h5` pattern (`n` is not a digit)
# 3. sort data by `n`
files = sorted(
    filter(lambda f: f.split(
        "/")[-1].split("_")[-1][:-3].isdigit(), glob.glob(os.path.join(data_dir, "*.h5"))),
    key=lambda filename: int(filename.split("/")[-1].split("_")[-1][:-3])
)

data = [h5_data(filename) for filename in files[::int(len(files)/10)]]

for i, dat in enumerate(data):
    plt.plot(dat.x, dat['u'], color=f"C{i}", label=f"iteration = {i}")

plt.legend()

plt.savefig(os.path.join(data_dir, "01.png"))

plt.title("Nagumo solution")
plt.show()
