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

name = "belousov_zhabotinsky_pirock"
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

    u_0 = mesh['fields']['u_0'][:]
    u_1 = mesh['fields']['u_1'][:]
    u_2 = mesh['fields']['u_2'][:]
    level = mesh['fields']['level'][:]

    index = np.argsort(centers)

    x = centers[index]
    u_0 = u_0[index]
    u_1 = u_1[index]
    u_2 = u_2[index]
    level = level[index]

    return x, (u_0, u_1, u_2), level, frame


data = sorted([read_file(file) for file in glob.glob(
    os.path.join(data_dir, "u_ite_*.h5"))], key=lambda tmp: tmp[-1])


def a_eq(b, c):
    f = 3.0
    q = 2e-4
    return (f*c)/(q+b)


fig, (ax_a, ax_b, ax_c, ax_level) = plt.subplots(4)

t_ini, t_end = 0., 1.
dt = (t_end - t_ini)/2000.

for i, idx in enumerate(range(0, len(data), int(len(data)/9))):
    a = data[idx][1][0]
    b = data[idx][1][1]
    c = data[idx][1][2]
    x = data[idx][0]

    ax_a.plot(x, a, label=f"$t^n = {idx*dt:.3f}$", color=f"C{i}")
    ax_a.plot(x, a_eq(b, c), "--", color=f"C{i}")
    ax_b.plot(x, b, color=f"C{i}")
    ax_c.plot(x, c, color=f"C{i}")

    ax_level.plot(data[idx][0], data[idx][2], color=f"C{i}")

ax_a.get_xaxis().set_ticklabels([])
ax_b.get_xaxis().set_ticklabels([])
ax_c.get_xaxis().set_ticklabels([])

ax_a.set_ylabel("a")
ax_b.set_ylabel("b")
ax_c.set_ylabel("c")
ax_level.set_ylabel("level")

ax_a.set_ylim(-50, 900)
ax_b.set_ylim(-1e-2, 0.9)
ax_c.set_ylim(-1e-3, 0.12)

ax_a.legend(loc="lower left", bbox_to_anchor=(
    0.01, 1.1), ncol=5, fontsize='x-small')

plt.savefig(os.path.join(img_dir, "01.png"), dpi=200)

if not arguments.only_save:
    plt.show()


if not arguments.only_save:
    fig, (ax_a, ax_b, ax_c, ax_level) = plt.subplots(4)

    a = data[0][1][0]
    b = data[0][1][1]
    c = data[0][1][2]
    x = data[0][0]
    sol_a = ax_a.plot(x, a, label="a", color="C0")[0]
    sol_a_eq = ax_a.plot(x, a_eq(b, c), "--", alpha=0.5,
                         label="a equilibre", color="C0")[0]
    sol_b = ax_b.plot(x, b, label="b", color="C1")[0]
    sol_c = ax_c.plot(x, c, label="c", color="C2")[0]

    ax_a.get_xaxis().set_ticklabels([])
    ax_b.get_xaxis().set_ticklabels([])
    ax_c.get_xaxis().set_ticklabels([])

    ax_a.set_ylabel("a")
    ax_b.set_ylabel("b")
    ax_c.set_ylabel("c")
    ax_level.set_ylabel("level")

    ax_a.set_ylim(-50, 850)
    ax_b.set_ylim(-1e-2, 0.8)
    ax_c.set_ylim(-1e-3, 0.12)

    levels = ax_level.plot(data[0][0], data[0][2])[0]

    def update(frame):
        a = data[frame][1][0]
        b = data[frame][1][1]
        c = data[frame][1][2]
        x = data[frame][0]

        sol_a.set_xdata(x)
        sol_a.set_ydata(a)

        sol_a_eq.set_xdata(x)
        sol_a_eq.set_ydata(a_eq(b, c))

        sol_b.set_xdata(x)
        sol_b.set_ydata(b)

        sol_c.set_xdata(x)
        sol_c.set_ydata(c)

        levels.set_xdata(x)
        levels.set_ydata(data[frame][2])

        return (sol_a, sol_b, sol_c, levels)

    N_frames = len(data)

    ani = animation.FuncAnimation(
        fig=fig, func=update, frames=N_frames, interval=1)

    plt.show()
