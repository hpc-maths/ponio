# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import glob
import argparse

import matplotlib.pyplot as plt
import numpy as np
import h5py

name = "nagumo_pirock"
data_dir = name+"_data"
img_dir = data_dir

parser = argparse.ArgumentParser(
    description=f"Compile, launch and plot results of example `{name}`")
parser.add_argument('--only-save', action='store_true',
                    help="Just save output in img directory")

arguments = parser.parse_args()

if arguments.only_save:
    img_dir = os.path.join("img", name)
    os.makedirs(img_dir, exist_ok=True)


def read_file(filename: str):
    # filename = os.path.join(data_dir, f"sol_1d_ite_{frame}.h5")
    frame = int(os.path.basename(filename).split("_")[-1][:-3])

    mesh = h5py.File(filename, 'r')['mesh']
    points = mesh['points']
    connectivity = mesh['connectivity']

    segments = np.zeros((connectivity.shape[0], 2, 2))
    segments[:, :, 0] = points[:][connectivity[:]][:, :, 0]

    centers = .5*(segments[:, 0, 0] + segments[:, 1, 0])
    index = np.argsort(centers)

    fields = {
        key: mesh['fields'][key][:][index]
        for key in mesh['fields'].keys()
    }

    x = centers[index]

    return x, fields, frame


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

data = sorted([read_file(file) for file in glob.glob(
    os.path.join(data_dir, "u_ite_*.h5"))], key=lambda tmp: tmp[-1])


fig, (ax_sol, ax_lvl) = plt.subplots(2)

for x, dat, i in data[::len(data)//9-1]:
    ax_sol.plot(x, dat['u'], color=f"C{i}", label=f"iteration = {i}")
    ax_lvl.plot(x, dat['level'], color=f"C{i}", label=f"iteration = {i}")

ax_sol.set_ylabel("solution")
ax_sol.legend(loc="upper left", fontsize='x-small')
ax_lvl.set_ylabel("levels")

plt.savefig(os.path.join(img_dir, "01.png"))

if not arguments.only_save:
    fig.suptitle("Nagumo solution")
    plt.show()
