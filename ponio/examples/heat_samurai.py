# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import glob
import argparse

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import h5py

name = "heat_samurai"
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


make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()


def read_file(filename: str):
    # filename = os.path.join(data_dir, f"sol_1d_ite_{frame}.h5")
    frame = int(os.path.basename(filename).split("_")[-1][:-3])

    mesh = h5py.File(filename, 'r')['mesh']
    points = mesh['points']
    connectivity = mesh['connectivity']

    segments = np.zeros((connectivity.shape[0], 2, 2))
    segments[:, :, 0] = points[:][connectivity[:]][:, :, 0]

    centers = .5*(segments[:, 0, 0] + segments[:, 1, 0])

    u = mesh['fields']['u'][:]
    level = mesh['fields']['level'][:]

    index = np.argsort(centers)

    x = centers[index]
    u = u[index]
    level = level[index]

    return x, u, level, frame


# plot solution and levels
plt.rcParams["figure.figsize"] = (8, 4)

data = sorted([read_file(file) for file in glob.glob(
    os.path.join(data_dir, "sol_1d_ite_*.h5"))], key=lambda tmp: tmp[-1])

fig, (ax1, ax2) = plt.subplots(2)

ax1.plot(data[0][0], data[0][1], color="C0", label="initial condition")
ax2.plot(data[0][0], data[0][2], color="C0", label="initial condition")

ax1.plot(data[-1][0], data[-1][1], color="C1", label="final time solution")
ax2.plot(data[-1][0], data[-1][2], color="C1", label="final time solution")

ax1.set_ylabel("solution")
ax2.set_ylabel("level")

ax1.legend()

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    fig.suptitle("Heat equation")
    plt.show()

# animation
fig, (ax1, ax2) = plt.subplots(2)

solution = ax1.plot(data[0][0], data[0][1])[0]
levels = ax2.plot(data[0][0], data[0][2])[0]


def update(frame):
    solution.set_xdata(data[frame][0])
    solution.set_ydata(data[frame][1])

    levels.set_xdata(data[frame][0])
    levels.set_ydata(data[frame][2])

    return solution


N_frames = len(data)

ani = animation.FuncAnimation(
    fig=fig, func=update, frames=N_frames, interval=10)

ani.save(filename=os.path.join(img_dir, "01.gif"), dpi=100, writer="pillow")

if not arguments.only_save:
    plt.show()
